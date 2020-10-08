import os
import abc
from collections import defaultdict

from Bio import SeqIO
from django.utils.functional import classproperty
import numpy
import regex

import azimuth3.model_comparison
import CFD

from ..helpers import reverse_complement, gc, extract_sequences, run_bowtie_multi, complement, dgn_to_regex
from . import Nuclease, OligoDict


class Cas9(Nuclease, abc.ABC):
    """Base Cas9 class. Can be subclassed for specific Cas9 variants."""

    @classproperty
    def _cut_site_position(cls):
        return cls.spacer_length - 3

    @classproperty
    def target_motif(cls):
        return f'(?P<upstream>.{{{cls.upstream_length}}})(?P<spacer>.{{{cls.spacer_length}}})(?P<PAM>{dgn_to_regex(cls.pam_motif)})(?P<downstream>.{{{cls.downstream_length}}})'

    @classmethod
    @abc.abstractmethod
    def _filter_offtarget(cls, seq):
        """Should return True, if seq is a valid off target"""

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
    def make_scaffold_oligos(cls, scaffold_name: str = None) -> OligoDict:

        scaffold = cls._scaffold_name_to_sequence(scaffold_name)
        return {
            'top': scaffold[5:-4],
            'bottom': reverse_complement(scaffold[9:])
        }

    @classmethod
    def make_spacer_oligos(cls, spacer_sequence, scaffold_name=None):
        target, scaffold = cls._spacer_to_cloning(spacer_sequence, scaffold_name)
        return {
            'top': 'cacc' + target + scaffold[:5].lower(),
            'bottom': reverse_complement(target + scaffold[:9].lower())
        }

    @classmethod
    def make_extension_oligos(cls, extension_sequence, scaffold_name=None):
        scaffold = cls._scaffold_name_to_sequence(scaffold_name)
        return {
            'top': scaffold[-4:].lower() + extension_sequence,
            'bottom': 'aaaa' + reverse_complement(extension_sequence)
        }

    @classmethod
    def make_nicking_oligos(cls, spacer_sequence, scaffold_name=None):
        target, scaffold = cls._spacer_to_cloning(spacer_sequence, scaffold_name)
        return {
            'top': 'cacc' + target,
            'bottom': reverse_complement(target + scaffold[:4].lower())
        }

    @classmethod
    def find_spacers(cls, reference_sequence, altered_sequence, start, end, spacer_search_range, **options):
        """Find candidate spacers for pegRNA selection.

        Finds all spacers with a cut site within spacer_search_range of the edit.
        Sorts spacers according to pam disruption, distance to edit and score.
        """
        spacers = []
        scoring_spacers = []
        sense = reference_sequence[:start + cls.downstream_from_cut_site]
        sense_offset = max(0, start - spacer_search_range - cls.cut_site_position)
        nucleotide_difference = len(altered_sequence) - len(reference_sequence)

        antisense = reverse_complement(
            reference_sequence[end - cls.downstream_from_cut_site - max(0, nucleotide_difference):])
        antisense_offset = end - max(0, nucleotide_difference)

        pam_motif = dgn_to_regex(cls.pam_motif) + '$'
        for match in regex.finditer(cls.target_motif,
                                    sense[-spacer_search_range - cls.cut_site_position - cls.downstream_from_cut_site:],
                                    regex.IGNORECASE, overlapped=True):
            spacer = match.group('spacer')
            pam = match.group('PAM')
            pos = sense_offset + match.start() + len(match.group('upstream'))
            pam_disrupted = not regex.search(pam_motif,
                                             altered_sequence[pos + len(spacer):pos + len(spacer) + len(pam)],
                                             regex.IGNORECASE)
            cut_site = pos + cls.cut_site_position - len(match.group('upstream'))
            distance = start - cut_site

            spacers.append({'spacer': spacer,
                            'position': pos,
                            'cut_site': cut_site,
                            'strand': 1,
                            'pam': (pam, pos + len(spacer)),
                            'pam_disrupted': pam_disrupted,
                            'distance': distance,
                            })

            scoring_spacers.append(match.group().upper())

        for match in regex.finditer(cls.target_motif, antisense[
                                                      -spacer_search_range - cls.cut_site_position - cls.downstream_from_cut_site:],
                                    regex.IGNORECASE, overlapped=True):
            spacer = match.group('spacer')
            pam = match.group('PAM')
            pos = antisense_offset + spacer_search_range - match.start() - len(
                match.group('upstream')) - 1 + cls.cut_site_position
            pam_disrupted = not regex.search(pam_motif, reverse_complement(altered_sequence[pos - len(spacer) - len(
                pam) + 1 + nucleotide_difference:pos + 1 + nucleotide_difference - len(spacer)]), regex.IGNORECASE)
            cut_site = pos - cls.cut_site_position + len(match.group('upstream')) + 1
            distance = cut_site - end + max(0, nucleotide_difference)

            spacers.append({'spacer': spacer,
                            'position': pos,
                            'cut_site': cut_site,
                            'strand': -1,
                            'pam': (pam, pos - len(spacer) - len(pam) + 1 + nucleotide_difference),
                            'pam_disrupted': pam_disrupted,
                            'distance': distance,
                            })

            scoring_spacers.append(match.group().upper())

        for i, score in enumerate(cls.score_spacers(scoring_spacers)):
            spacers[i]['score'] = score

        return sorted(spacers, key=lambda x: (not x['pam_disrupted'], x['distance'], x['score']))

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
        bowtie_out = os.path.join(write_folder, 'bt_out')

        with open(spacers_file, 'w') as f:
            for spacer in spacer_sequences:
                f.write(f'{spacer}\n')

        # Run bowtie on the spacer sequences and parse output for extracting spacer and PAM from 2bit file
        run_bowtie_multi(spacers_file, assembly, bowtie_out)
        parsed = False
        with open(bowtie_out) as bt:
            f = open(in_file, 'w')
            for i, line in enumerate(bt):
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

                if (i+1) % 100000 == 0:
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
                                    {'off target site': match_to['mismatched_seq'][record.name], 'pam': seq[-len(cls.pam_motif):].upper(),
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

    @classmethod
    def find_nicking_spacers(cls, reference_sequence, altered_sequence, spacer_strand, cut_site, scaffold,
                             nicking_range, **options):
        """Find spacers for nicking the opposite strand."""
        spacers = []
        scoring_spacers = []

        nt_difference = len(altered_sequence) - len(reference_sequence)

        if spacer_strand == 1:
            reference_sequence = reverse_complement(reference_sequence)
            altered_sequence = reverse_complement(altered_sequence)
            cut_site = len(altered_sequence) - cut_site
        sequence = altered_sequence[
                   max(0, cut_site - cls.cut_site_position - nicking_range):
                   min(len(altered_sequence), cut_site + cls.downstream_from_cut_site + nicking_range)
                   ].upper()
        ref = reference_sequence[
              max(0, cut_site - cls.cut_site_position - nicking_range):
              min(len(reference_sequence), cut_site + cls.downstream_from_cut_site + nicking_range - nt_difference)
              ].upper()

        cut_site = cls.cut_site_position + nicking_range

        for match in regex.finditer(cls.target_motif, sequence, regex.IGNORECASE, overlapped=True):
            spacer = match.group('spacer')
            pos = match.start() + len(match.group('upstream'))
            wt_pos = pos
            cut = match.start() + cls.cut_site_position
            nick_location = cut_site - cut
            if nick_location < 0:
                wt_pos -= nt_difference

            kind = '3'
            wt_score = 1
            alt_bind = sequence[pos:pos + len(spacer) + len(cls.pam_motif)].upper()
            wt_bind = ref[wt_pos:wt_pos + len(spacer) + len(cls.pam_motif)].upper()
            if cls._is3b(alt_bind, wt_bind):
                kind = '3b'
                wt_score = cls._calc_wt_score(alt_bind, wt_bind)

            info = cls.make_nicking_oligos(spacer, scaffold)
            info['position'] = nick_location
            info['spacer'] = spacer
            info['kind'] = kind
            info['wt_score'] = wt_score
            info['offset'] = cls.cut_site_position - len(match.group('upstream'))

            spacers.append(info)
            scoring_spacers.append(match.group().upper())

        for i, score in enumerate(cls.score_spacers(scoring_spacers)):
            spacers[i]['score'] = score

        return sorted(spacers, key=lambda x: (x['wt_score'], not (abs(x['position']) > 50), -x['score']))

    @classmethod
    def _make_pbs_sequence(cls, reference, pbs_min_length, pbs_max_length, **options):
        """Find a suggested PBS length, and generate all possible PBS candidate lengths.

        Selects the shortest PBS sequence with a GC content in the range [0.4,0.6].
        If no sequence is within this range, selects the shortest PBS with a GC content closest to 0.5.
        """
        pbs_length = pbs_min_length - 1
        lengths = []
        while pbs_length < pbs_max_length:
            pbs_length += 1
            pbs = reference[-pbs_length:]
            if 0.4 <= gc(pbs) <= 0.6:
                break
            lengths.append((abs(0.5 - gc(pbs)), len(pbs), pbs))
        else:
            pbs = sorted(lengths, key=lambda x: x[:1])[0][2]

        # Create all possible PBS sequences within range limits.
        alt_lengths = [reference[-pbs_length:] for pbs_length in range(pbs_min_length, pbs_max_length + 1)]
        return pbs, alt_lengths

    @classmethod
    def _make_rt_sequence(cls, reference, cut_dist, nucleotide_difference, alteration_length, rt_min_length,
                          rt_max_length, **options):
        """Find a suggested RT template length, and generate alterniative RT tempalte lengths."""
        rt_template_length = rt_min_length
        to_position = cut_dist + rt_template_length + nucleotide_difference
        rt_template = reference[:to_position]
        # Extension sequence should not start with a 'C'
        # Increase reverse transcription template until it doesn't end with a 'G'
        # For large alterations, longer template is probably preferred

        while rt_template_length <= alteration_length * 2 and rt_template_length <= rt_max_length:
            rt_template += reference[to_position]
            to_position += 1
            rt_template_length += 1

        # Create all possible RT templates within range limits that do not end with a 'G'.
        lengths = []
        for rt_template_length in range(rt_min_length, rt_max_length + 1):
            template = reference[:cut_dist + rt_template_length + alteration_length]
            lengths.append(template)
        return rt_template, lengths

    @classmethod
    def make_extension_sequence(cls, reference_sequence, altered_sequence, spacer_strand, spacer_cut_site, cut_dist,
                                alteration_length, pbs_min_length, pbs_max_length, rt_min_length, rt_max_length,
                                **options):
        """Create the pegRNA extension sequence.

        pegRNA extension sequences consist of a PBS and a RT template.
        The PBS is upstream of the cut site. The RT template is downstream of the cut site and contains the edit sequence.

        """
        nucleotide_difference = len(altered_sequence) - len(reference_sequence)
        if spacer_strand == 1:
            nucleotide_difference = min(0, nucleotide_difference)
            pbs_reference = reference_sequence[:spacer_cut_site]
            rt_reference = altered_sequence[spacer_cut_site:]
        else:
            pbs_reference = reverse_complement(reference_sequence[spacer_cut_site:])
            rt_reference = reverse_complement(altered_sequence[:spacer_cut_site + nucleotide_difference])
        pbs, pbs_lengths = cls._make_pbs_sequence(pbs_reference.upper(), pbs_min_length, pbs_max_length)
        rt, rt_lengths = cls._make_rt_sequence(rt_reference, cut_dist, nucleotide_difference, alteration_length,
                                               rt_min_length, rt_max_length)

        pbs_length = len(pbs)
        rt_length = len(rt)

        # Generate all combinations of PBS and RT template sequences that are not identical to the primary suggestion.
        alternate_extensions = []
        for alt_pbs in pbs_lengths:
            alt_pbs_length = len(alt_pbs)
            alt_pbs_gc = round(gc(alt_pbs), 2)
            for alt_rt in rt_lengths:
                alt_rt_length = len(alt_rt)
                alt_rt_gc = round(gc(alt_rt), 2)
                if alt_pbs_length == pbs_length and alt_rt_length == rt_length:
                    continue
                alternate_extensions.append({'pbs_length': alt_pbs_length, 'rt_template_length': alt_rt_length,
                                             'sequence': reverse_complement(alt_pbs + alt_rt), 'pbs_gc': alt_pbs_gc,
                                             'rt_gc': alt_rt_gc})

        return pbs_length, rt_length, reverse_complement(pbs + rt), alternate_extensions


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
    def _make_rt_sequence(cls, reference, cut_dist, nucleotide_difference, alteration_length, rt_min_length,
                          rt_max_length, **options):
        """Find a suggested RT template length, and generate alterniative RT tempalte lengths."""
        rt_template_length = rt_min_length
        to_position = cut_dist + rt_template_length + nucleotide_difference
        rt_template = reference[:to_position]
        last_valid = rt_template
        # Extension sequence should not start with a 'C'
        # Increase reverse transcription template until it doesn't end with a 'G'
        # For large alterations, longer template is probably preferred

        while (rt_template.endswith(
                'G') or rt_template_length <= alteration_length * 2) and rt_template_length <= rt_max_length:
            rt_template += reference[to_position]
            if not rt_template.endswith('G'):
                last_valid = rt_template
            to_position += 1
            rt_template_length += 1

        # If canceled by cutoff, might end with G .. Take the last valid length.
        rt_template = last_valid

        # Create all possible RT templates within range limits that do not end with a 'G'.
        lengths = []
        for rt_template_length in range(rt_min_length, rt_max_length + 1):
            template = reference[:cut_dist + rt_template_length + alteration_length]
            if not template.endswith('G'):
                lengths.append(template)
        return rt_template, lengths


class SpCas9(SpCas9Base):
    """SpCas9 nuclease."""

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
    """SpCas9-NG"""

    pam_motif = 'NGN'

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq[-2] != 'G'


class SaCas9(Cas9):
    """SaCas9"""

    pam_motif = 'NNGRRT'

    spacer_length = 21

    scaffolds = {
        'sgRNA':    'GTTTTAGTACTCTGGAAACAGAATCTACTAAAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGA',
        'sgRNAmod': 'GTTCTAGTACTCTGGAAACAGAATCTACTAGAACAAGGCAAAATGCCGTGTTTATCTCGTCAACTTGTTGGCGAGA',
    }
    default_scaffold = 'sgRNA'

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq[-4:-1] in ['GAA', 'GAG', 'GGA', 'GGG']


class SaCas9KKH(SaCas9):
    """SaCas9KKH"""

    pam_motif = 'NNNRRT'

    @classmethod
    def _filter_offtarget(cls, seq):
        return seq[-3:-1] in ['AA', 'AG', 'GA', 'GG']


class CjCas9(Cas9):
    """CjCas9"""

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
