import os
import abc
import subprocess
import distutils.util

import django
from django.utils.functional import classproperty
import numpy
import regex

import azimuth3.model_comparison
import CFD

from design.helpers import complement, dgn_to_regex
from .. import Nuclease
from .plasmids import GGAssembly, LibraryCloning, Synthetic
from .design_rules import Anzalone


class Cas9(Nuclease, abc.ABC):
    """Base Cas9 class. Can be subclassed for specific Cas9 variants."""

    cloning_strategies = {
        'pegRNA-GG-acceptor': GGAssembly,
        'Library': LibraryCloning,
        'Synthetic': Synthetic,
    }

    design_strategies = {
        'Anzalone': Anzalone
    }

    @classproperty
    def _cut_site_position(cls):
        return cls.spacer_length - 3

    @classproperty
    def target_motif(cls):
        return f'(?P<upstream>[ACTG]{{{cls.upstream_length}}})(?P<spacer>[ACTG]{{{cls.spacer_length}}})(?P<PAM>{dgn_to_regex(cls.pam_motif)})(?P<downstream>[ACTG]{{{cls.downstream_length}}})'

    @classmethod
    def _is3b(cls, alt_bind, wt_bind):
        """Should return True if the sgRNA only binds to the altered sequence"""
        if regex.match(fr'{alt_bind[:cls.spacer_length]}{dgn_to_regex(cls.pam_motif)}', wt_bind):
            return False
        return True

    @classmethod
    def _calc_wt_score(cls, alt_bind, wt_bind):
        return 1

    @classmethod
    def make_full_pegrna(cls, spacer, scaffold, extension):
        return ''.join([spacer, scaffold, extension])

    @classmethod
    def find_off_targets(cls, spacer_sequences, assembly, write_folder):
        """Find off targets for a list of spacer sequences."""
        mm_pattern = regex.compile(r'(?P<position>\d+):(?P<to>\w)>(?P<from>\w)')
        counts = []
        binders = []

        for _ in spacer_sequences:
            counts.append([0] * 4)
            binders.append([])

        spacers_file = os.path.join(write_folder, 'spacers')

        with open(spacers_file, 'w') as f:
            for spacer in spacer_sequences:
                f.write(f'{spacer}{"N"*len(cls.pam_motif)}\n')

        # Run bowtie on the spacer sequences and parse output for spacer and PAM sequence
        cmd = f"{django.conf.settings.DESIGN_BOWTIE_PATH} -n 3 -l {cls.spacer_length} -e {(3 + len(cls.pam_motif))*30} -a --best --sam-nohead -x {django.conf.settings.DESIGN_BOWTIE_GENOMES_FOLDER}/{assembly}/{assembly} --suppress 5,6,7 -r {spacers_file} -y --threads {django.conf.settings.DESIGN_BOWTIE_THREADS}"
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
        for i, line in enumerate(iter(p.stdout.readline, '')):
            line = line.strip('\n')
            index, strand, scaffold, position, mismatches = line.split('\t')
            position = int(position)
            index = int(index)
            mismatches = regex.findall(mm_pattern, mismatches)
            match = list(spacer_sequences[index].upper()) + ['N'] * len(cls.pam_motif)
            if strand == '+':
                for mm in mismatches:
                    match[int(mm[0])] = mm[1].lower()
            else:
                position += len(cls.pam_motif)
                for mm in mismatches:
                    match[int(mm[0])] = complement(mm[1].lower())
            match = ''.join(match)

            spacer, pam = match[:-len(cls.pam_motif)], match[-len(cls.pam_motif):]
            pam = pam.upper()

            if cls._filter_offtarget(pam):
                mm_count = len(mismatches) - len(cls.pam_motif)
                counts[index][mm_count] += 1
                if len(binders[index]) < 20:
                    binders[index].append({
                        'off target site': spacer,
                        'pam': pam,
                        'chr': scaffold,
                        'position': position,
                        'strand': strand,
                        'mismatches': mm_count})

        return counts, binders


class SpCas9Base(Cas9, abc.ABC):
    """Base class for SpCas9 variants"""

    scaffolds = {
        'sgRNA': 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC',
        'sgRNA2.0': 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC',
        'sgRNA2.1': 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTTGCTGGAAACAGCAAAGTGGCACCGAGTCGGTGC',
    }

    default_scaffold = 'sgRNA'
    spacer_length = 20

    _options = {
        'remove extensions starting with C': (bool, True, None, 'SpCas9 pegRNA extensions that start with a C may be less efficient'),  # Filter out extensions that start with a C
    }

    @classmethod
    def filter_extension(cls, seq, **options):
        try:
            remove_C = bool(distutils.util.strtobool(str(options['remove extensions starting with C'])))
        except KeyError:
            remove_C = cls._options['remove extensions starting with C'][1]
        return remove_C and seq.endswith('G')


class SpCas9(SpCas9Base):
    """SpCas9 nuclease.

    20nt spacer with NGG PAM.

    Filters out RT-templates that end with a 'G'

    On-target scores based on Doench. et al, 2016.

    Nicking sgRNAs scored using CFD score against wild type sequence, Doench. et al, 2014.
    """

    pam_motif = 'NGG'

    upstream_length = 4
    downstream_length = 3
    can_score_spacers = True

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq[-2:] in ['AG', 'GG']

    @classmethod
    def _calc_wt_score(cls, alt_bind, wt_bind):
        return CFD.calc_cfd(alt_bind, wt_bind[:-3], wt_bind[-2:])

    @classmethod
    def score_spacers(cls, spacers):
        """Use azimuth to score spacers."""
        if spacers:
            return azimuth3.model_comparison.predict(numpy.array(spacers), None, None, silent=True)
        return []


class SpCas9NG(SpCas9Base):
    """SpCas9-NG

    20 nt spacer with NGN PAM.

    Filters out RT-templates that end with a 'G'
    """

    pam_motif = 'NGN'

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq[-2] == 'G'


class SaCas9(Cas9):
    """SaCas9

    21nt spacer with NNGRRT PAM.
    """

    pam_motif = 'NNGRRT'

    spacer_length = 21

    scaffolds = {
        'sgRNA': 'GTTTTAGTACTCTGGAAACAGAATCTACTAAAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGA',
        'sgRNAmod': 'GTTCTAGTACTCTGGAAACAGAATCTACTAGAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGA',
    }
    default_scaffold = 'sgRNA'

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq[-4:-1] in ['GAA', 'GAG', 'GGA', 'GGG']


class SaCas9KKH(SaCas9):
    """SaCas9KKH

    21nt spacer with NNNRRT PAM.
    """

    pam_motif = 'NNNRRT'

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq[-3:-1] in ['AA', 'AG', 'GA', 'GG']


class CjCas9(Cas9):
    """CjCas9

    22nt spacer with NNNNACAC PAM.
    """

    pam_motif = 'NNNNACAC'

    spacer_length = 22

    scaffolds = {
        'sgRNA': 'GTTTTAGTCCCTGAAGGGACTAAAATAAAGAGTTTGCGGGACTCTGCGGGGTTACAATCCCCTAAAACCGC',
        'sgRNAmod': 'GTTCTAGTCCCTGAAGGGACTAGAATAAAGAGTTTGCGGGACTCTGCGGGGTTACAATCCCCTAAAACCGC',
    }

    default_scaffold = 'sgRNA'

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq.endswith('ACAC')
