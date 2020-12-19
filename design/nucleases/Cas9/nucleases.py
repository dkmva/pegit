import os
import abc
from collections import defaultdict

from Bio import SeqIO
from django.utils.functional import classproperty
import numpy
import regex

import azimuth3.model_comparison
import CFD

from design.helpers import reverse_complement, extract_sequences, run_bowtie_multi, complement, dgn_to_regex
from .. import Nuclease
from .plasmids import GGAssembly, LibraryCloning
from .design_rules import Anzalone


class Cas9(Nuclease, abc.ABC):
    """Base Cas9 class. Can be subclassed for specific Cas9 variants."""

    cloning_strategies = {
        'pegRNA-GG-acceptor': GGAssembly,
        'Library': LibraryCloning,
    }

    design_strategies = {
        'Anzalone': Anzalone
    }

    @classmethod
    def do_alternate_cloning(cls, strategy=None, scaffold=None, **options):
        if strategy is None:
            strategy = next(iter(cls.cloning_strategies.values()))
        if isinstance(strategy, str):
            strategy = cls.cloning_strategies[strategy]
        scaffold = cls._scaffold_name_to_sequence(scaffold)
        return strategy.alternate_extension(scaffold=scaffold, **options)

    @classproperty
    def _cut_site_position(cls):
        return cls.spacer_length - 3

    @classproperty
    def target_motif(cls):
        return f'(?P<upstream>.{{{cls.upstream_length}}})(?P<spacer>.{{{cls.spacer_length}}})(?P<PAM>{dgn_to_regex(cls.pam_motif)})(?P<downstream>.{{{cls.downstream_length}}})'

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
    def find_off_targets(cls, spacer_sequences, assembly, write_folder):
        """Find off targets for a list of spacer sequences."""
        mm_pattern = regex.compile(r'(?P<position>\d+):(?P<to>\w)>(?P<from>\w)')
        spacer_len = len(spacer_sequences[0])
        counts = []
        binders = []

        for _ in spacer_sequences:
            counts.append([0] * 4)
            binders.append([])

        match_to = {
            'strand': {}, 'mismatches': defaultdict(dict), 'mismatched_seq': {}, 'position': {}, 'scaffold': {}
        }

        in_file = os.path.join(write_folder, 'extract_sequences_in')
        out_file = os.path.join(write_folder, 'extract_sequences_out')
        spacers_file = os.path.join(write_folder, 'spacers')

        with open(spacers_file, 'w') as f:
            for spacer in spacer_sequences:
                f.write(f'{spacer}\n')

        # Run bowtie on the spacer sequences and parse output for extracting spacer and PAM from 2bit file
        bowtie_out = run_bowtie_multi(spacers_file, assembly)
        parsed = False
        # with open(bowtie_out) as bt:
        f = open(in_file, 'w')
        for i, line in enumerate(bowtie_out):
            parsed = False
            line = line.strip('\n')
            index, strand, scaffold, position, mismatches = line.split('\t')
            position = int(position)
            index = int(index)
            mismatches = regex.findall(mm_pattern, mismatches)

            match = list(spacer_sequences[index].upper())
            if strand == '+':
                position_string = f'{scaffold}:{position}-{position + spacer_len + len(cls.pam_motif)}'
                for mm in mismatches:
                    match[int(mm[0])] = mm[1].lower()
            else:
                if (position - len(cls.pam_motif)) < 0:
                    continue
                position_string = f'{scaffold}:{position - len(cls.pam_motif)}-{position + spacer_len}'
                for mm in mismatches:
                    match[int(mm[0])] = complement(mm[1].lower())
            f.write(f'{position_string}\n')
            match = ''.join(match)

            match_to['strand'][position_string] = strand
            match_to['mismatches'][position_string][index] = len(mismatches)
            match_to['position'][position_string] = position
            match_to['scaffold'][position_string] = scaffold
            match_to['mismatched_seq'][position_string] = match

            if (i + 1) % 100000 == 0:
                f.close()
                # Get the sequences from 2bit file, filter out hits without a potential PAM
                extracted_sequences = SeqIO.parse(extract_sequences(in_file, assembly, out_file), 'fasta')
                for record in extracted_sequences:
                    strand = match_to['strand'][record.name]
                    seq = str(record.seq)
                    if strand == '-':
                        seq = reverse_complement(seq)
                    if not cls._filter_offtarget(seq):
                        continue
                    # Only first 20 hits are returned, all hits are counted.
                    for idx in match_to['mismatches'][record.name].keys():
                        counts[idx][match_to['mismatches'][record.name][idx]] += 1
                        if len(binders[idx]) < 20:
                            binders[idx].append(
                                {'off target site': match_to['mismatched_seq'][record.name],
                                 'pam': seq[-len(cls.pam_motif):].upper(),
                                 'chr': match_to['scaffold'][record.name],
                                 'position': match_to['position'][record.name],
                                 'strand': strand, 'mismatches': match_to['mismatches'][record.name][idx]})
                match_to = {
                    'strand': {}, 'mismatches': defaultdict(dict), 'mismatched_seq': {}, 'position': {},
                    'scaffold': {}
                }
                parsed = True
                f = open(in_file, 'w')

        f.close()
        if not parsed:
            extracted_sequences = SeqIO.parse(extract_sequences(in_file, assembly, out_file), 'fasta')
            for record in extracted_sequences:
                strand = match_to['strand'][record.name]
                seq = str(record.seq)
                if strand == '-':
                    seq = reverse_complement(seq)
                if not cls._filter_offtarget(seq):
                    continue
                # Only first 20 hits are returned, all hits are counted.
                for idx in match_to['mismatches'][record.name].keys():
                    counts[idx][match_to['mismatches'][record.name][idx]] += 1
                    if len(binders[idx]) < 20:
                        binders[idx].append(
                            {'off target site': match_to['mismatched_seq'][record.name],
                             'pam': seq[-len(cls.pam_motif):].upper(),
                             'chr': match_to['scaffold'][record.name],
                             'position': match_to['position'][record.name],
                             'strand': strand, 'mismatches': match_to['mismatches'][record.name][idx]})

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

    @classmethod
    def filter_extension(cls, seq):
        return seq.endswith('G')


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
        return seq[-2] != 'G'


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
