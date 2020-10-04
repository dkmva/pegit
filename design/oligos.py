"""
Module for designing cloning oligos.
"""
from collections import defaultdict
import copy
import regex

import django

from design.helpers import is_valid_dna, is_valid_degenerate, reverse_complement, degenerate_to_nucleotides, \
    run_primer3
from design.nucleases import NUCLEASES


class OligoSet:
    """Set of oligos for cloning pegRNAs."""

    def __init__(self, tracker, spacer, scaffold=None, nuclease=None, repair=False):
        if nuclease is None:
            nuclease = django.conf.settings.DESIGN_CONF['default_nuclease']
        if isinstance(nuclease, str):
            nuclease = NUCLEASES[nuclease]
        if scaffold is None:
            self.scaffold = nuclease.default_scaffold

        self.spacer_sequence = spacer['spacer']
        self.spacer_position = spacer['position']
        self.spacer_cut_site = spacer['cut_site']
        self.spacer_strand = spacer['strand']
        self.pam = spacer['pam']
        self.pam_disrupted = spacer['pam_disrupted']
        self.spacer_distance = spacer['distance']
        self.spacer_score = spacer['score']
        self.visual_spacer = spacer['visual_spacer']
        self.pam_silenced = False
        self.nuclease = nuclease
        self.tracker = tracker
        self.pbs_length = None
        self.rt_template_length = None
        self.oligos = {}
        self.nicking_spacers = {}
        self.alternate_extension = []
        self.repair = repair

    @property
    def info(self):
        info = {
            'spacer': self.spacer_sequence,
            'score': self.spacer_score,
            'position': self.spacer_position,
            'pam_disrupted': self.pam_disrupted,
            'pam_silenced': self.pam_silenced,
            'distance': self.spacer_distance,
            'strand': self.spacer_strand,
            'nuclease': self.nuclease.__name__,
            'pbs_length': self.pbs_length,
            'rt_template_length': self.rt_template_length,
            'oligos': self.oligos,
            'extension': self.extension,
            'visual_extension': self.visual_exension,
            'visual_spacer': self.visual_spacer,
            'nicking': self.nicking_spacers,
            'alternate_extensions': self.alternate_extension,
            'pbs': self.pbs,
            'rt_template': self.rt_template,
        }
        return info

    def make_oligos(self, degenerate_sequence, silence_pam=False, **options):
        """Make oligos for cloning pegRNAs"""
        if not self.pam_disrupted and silence_pam:
            self.tracker = self.tracker.copy()
            if self.spacer_strand == 1:
                pam = self.nuclease.pam_motif
            else:
                pam = reverse_complement(self.nuclease.pam_motif)

            for i, j in enumerate(range(self.pam[1], self.pam[1] + len(self.pam[0]))):
                dgn = degenerate_sequence[j]
                pam_nt = pam[i]
                pam_dgn = degenerate_to_nucleotides[pam_nt]
                nt_dgn = degenerate_to_nucleotides[dgn]
                if pam_dgn < nt_dgn:
                    for nt in nt_dgn:
                        if nt not in pam_dgn:
                            self.tracker.substitute(nt.lower(), j)
                            self.pam_silenced = True
                            break

        reference_sequence = self.tracker.original_sequence
        altered_sequence = str(self.tracker)
        self.oligos['spacer'] = self.nuclease.make_spacer_oligos(self.spacer_sequence, self.scaffold)
        self.oligos['scaffold'] = self.nuclease.make_scaffold_oligos(self.scaffold)
        self._make_extension_sequence(**options)
        self.oligos['extension'] = self.nuclease.make_extension_oligos(self.extension, self.scaffold)
        for extension in self.alternate_extension:
            extension['oligos'] = self.nuclease.make_extension_oligos(extension['sequence'], self.scaffold)

        if self.repair:
            reference_sequence, altered_sequence = altered_sequence, reference_sequence
        spacers = self.nuclease.find_nicking_spacers(
            reference_sequence,
            altered_sequence,
            self.spacer_strand,
            self.spacer_cut_site,
            self.scaffold, **options)

        for spacer in spacers:

            visual_spacer = spacer['spacer']
            position = spacer['position']
            if self.spacer_strand == 1:
                visual_spacer = reverse_complement(visual_spacer)

            pos = self.spacer_cut_site
            if self.spacer_strand == 1:
                pos += position - len(spacer['spacer']) + spacer['offset']
            else:
                pos -= position + spacer['offset']

            if self.repair:
                visual_spacer = self.tracker.seq_from_original_coordinates(pos, pos + len(visual_spacer))
                spacer['push'] = 0
                if self.tracker.number_of_insertions:
                    spacer['push'] = (pos - self.tracker.index[pos][0]) + (self.spacer_cut_site -
                                self.tracker.index[self.spacer_cut_site][0])
            else:
                visual_spacer = self.tracker.seq_from_new_coordinates(pos, pos + len(visual_spacer))
                spacer['push'] = 0
                if self.tracker.number_of_deletions:
                    spacer['push'] = (self.tracker.index[pos][0] - pos) + (self.tracker.index[self.spacer_cut_site][0] - self.spacer_cut_site)
            if self.spacer_strand == 1:
                spacer['push'] -= self.tracker.number_of_insertions + self.tracker.number_of_deletions

            if self.spacer_strand == 1:
                visual_spacer = reverse_complement(visual_spacer)

            spacer['visual_spacer'] = visual_spacer

        self.nicking_spacers = spacers
        return self.oligos

    def _make_extension_sequence(self, **options):
        """Make the extension sequence for the pegRNA.

        Also creates visual_extension, which can be used by the frontend for visualization.
        """

        reference_sequence = self.tracker.original_sequence
        altered_sequence = str(self.tracker)
        alteration_position = sum(self.tracker.alterations()[0])
        if alteration_position >= self.spacer_cut_site:
            cut_dist = alteration_position+self.tracker.alteration_length - self.spacer_cut_site
        else:
            cut_dist = self.spacer_cut_site - alteration_position

        if self.repair:
            reference_sequence, altered_sequence = altered_sequence, reference_sequence
        pbs_length, rt_template_length, extension, alternate_extensions = self.nuclease.make_extension_sequence(
            reference_sequence,
            altered_sequence,
            self.spacer_strand,
            self.spacer_cut_site,
            cut_dist,
            self.tracker.number_of_alterations,
            **options)

        self.pbs_length = pbs_length
        self.rt_template_length = rt_template_length
        self.extension = extension
        self.rt_template = extension[:rt_template_length]
        self.pbs = extension[rt_template_length:]
        self.alternate_extension = alternate_extensions
        visual_extension = extension

        start = self.spacer_cut_site - self.pbs_length
        if self.spacer_strand == 1:
            visual_extension = reverse_complement(visual_extension)

        if self.repair:
            if self.spacer_strand == -1:
                start = self.spacer_cut_site - self.rt_template_length + self.tracker.number_of_deletions - self.tracker.number_of_insertions
            visual_extension = self.tracker.seq_from_original_coordinates(start, start + len(visual_extension))
        else:
            if self.spacer_strand == -1:
                start = self.spacer_cut_site - self.rt_template_length - self.tracker.number_of_deletions + self.tracker.number_of_insertions
            visual_extension = self.tracker.seq_from_new_coordinates(start, start + len(visual_extension))
        if self.spacer_strand == 1:
            visual_extension = reverse_complement(visual_extension)

        self.visual_exension = visual_extension


class SequenceObject:
    """Basic sequence object."""
    def __init__(self, sequence, organism, name=None, id=None, source=None):
        self.sequence = sequence
        self.organism = organism
        self.upstream = ''
        self.downstream = ''
        self.name = name
        self.id = id
        self.source = source

    def contextualize_alteration(self, alteration):
        return alteration

    @property
    def degenerate_sequence(self):
        return []

    @property
    def codon_usage_table(self):
        return self.organism.codon_usage_table


class AlterationTracker:
    """Sequence alteration tracker.

    A object to keep track of modifications to a DNA sequence.

    """
    def __init__(self, sequence):
        #if not is_valid_dna(sequence):
        #    raise ValueError(f'Sequence is not valid DNA {sequence}')

        sequence = str(sequence).upper()
        self._original_sequence = sequence
        self.sequence_range = range(len(sequence))
        self.sequence_range = range(-1, len(sequence)+1)

        self.sequence = {i: [[nucleotide, nucleotide]] for i, nucleotide in enumerate(sequence)}
        self.sequence[-1] = []
        self.sequence[len(sequence)] = []
        self.sequence = dict(sorted(self.sequence.items()))
        self.index = [(i, 0) for i in range(len(sequence))]
        self.alteration_list = {i: ['n'] for i in range(len(sequence))}  # No longer a list ..
        self.alteration_list[-1] = []  # For insertions a first position
        self.alteration_list[len(sequence)] = []  # For insertions at last position
        self.alteration_list = dict(sorted(self.alteration_list.items()))
        self.degeneracy = []

    def __str__(self):
        sequence = []
        for pos, alts in zip(self.sequence.values(), self.alteration_list.values()):
            for nts, alt in zip(pos, alts):
                if alt == 'd':
                    continue
                if alt == 'n':
                    sequence.append(nts[0])
                else:
                    sequence.append(nts[0].lower())

        return ''.join(sequence)

    def __iter__(self):
        for i, j in self.index:
            yield self.sequence[i][j][0]

    def __getitem__(self, index):
        try:
            i, j = index
            return self.sequence[i][j][0]
        except TypeError:
            pass
        try:
            index = slice(index, index + 1)
        except TypeError:
            pass
        return ''.join([self.sequence[i][j][0] for i, j in self.index[index]])

    @property
    def original_sequence(self):
        """Return the original sequence, with altered positions as lowercase."""
        seq = []
        for i, nt in enumerate(self._original_sequence):
            if self.alteration_list[i][0] != 'n':
                seq.append(nt.lower())
            else:
                seq.append(nt)
        return ''.join(seq)

    def substitute(self, nucleotide, position, degenerate=None, in_place=True, force=False):
        """Make a nucleotide substitution.

                Required arguments:
                nucleotide -- new nucleotide at position
                position -- position to substitute

                Optional arguments:
                degenerate -- degeneracy allowed at the position (if any), defaults to only new nucleotide
                in_place -- if True, change the tracker instance, if False return a modified copy, defaults to True.
                force -- if True, forces DNA change, even if disallowed by degeneracy, defaults to False.
                """
        if not (len(nucleotide) == 1 and is_valid_dna(nucleotide)):
            raise ValueError(f'{nucleotide} is not a single nucleotide')
        nucleotide = nucleotide.upper()
        i, j = self.index[position]
        new = self if in_place else self.copy()

        if degenerate is None:
            degenerate = nucleotide

        for dgn in new.degeneracy:
            if i in dgn:
                if force:
                    dgn.substitute(position, degenerate)
                elif nucleotide.upper() not in dgn[position]:
                    return False

        new.sequence[i][j][0] = nucleotide
        if j == 0:
            if new.sequence[i][j][1].upper() == nucleotide.upper():
                new.alteration_list[i][j] = 'n'
            else:
                new.alteration_list[i][j] = 'm'
        new.remake_index()
        return new

    def insert(self, nucleotide, position, degenerate=None, in_place=True, force=False):
        """Insert a nucleotide.

                Required arguments:
                nucleotide -- new nucleotide at position
                position -- position of insertion

                Optional arguments:
                degenerate -- degeneracy allowed at the position (if any), defaults to only new nucleotide
                in_place -- if True, change the tracker instance, if False return a modified copy, defaults to True.
                force -- if True, forces DNA change, even if disallowed by degeneracy, defaults to False.
                """
        nucleotide = nucleotide.upper()
        i, j = self.index[position]
        new = self if in_place else self.copy()

        if degenerate is None:
            degenerate = nucleotide

        for dgn in new.degeneracy:
            if position in dgn:
                if force:
                    dgn.insert(position, degenerate)
                else:
                    return False

        if j > 0:
            new.sequence[i].insert(j, [nucleotide, '-'])
            new.alteration_list[i].insert(j, 'i')
        else:
            new.sequence[i - 1].append([nucleotide, '-'])
            new.alteration_list[i - 1].append('i')

        new.remake_index()
        return new

    def delete(self, position, j=None, in_place=True, force=False):
        """Delete a nucleotide.

                Required arguments:
                position -- position to delete

                Optional arguments:
                j -- can be used to specify position by tracker index.
                in_place -- if True, change the tracker instance, if False return a modified copy, defaults to True.
                force -- if True, forces DNA change, even if disallowed by degeneracy, defaults to False.
                """
        i = position
        if j is None:
            i, j = self.index[i]
        new = self if in_place else self.copy()

        for dgn in new.degeneracy:
            if position in dgn:
                if force:
                    dgn.delete(position)
                else:
                    return False

        if len(new.sequence[i]) > 1:
            del new.sequence[i][j]
            del new.alteration_list[i][j]
        else:
            new.sequence[i][j][0] = '-'
            new.alteration_list[i] = ['d']

        new.remake_index()
        return new

    def remake_index(self):
        self.index = [(i, j) for i in self.sequence_range
                      for j in range(len(self.sequence[i]))
                      if self.alteration_list[i][j] != "d"]

    def copy(self):
        return copy.deepcopy(self)

    def alteration_string(self):
        return ''.join([alteration for alterations in self.alteration_list.values() for alteration in alterations])

    def alterations(self):
        return [(i, j) for i in self.sequence_range
                for j in range(len(self.sequence[i]))
                if self.alteration_list[i][j] != 'n']

    def alignment_string(self):
        conversion = {'n': '|', 'm': '.', 'd': ' ', 'i': ' '}
        return ''.join([conversion[alteration] for alteration in self.alteration_string()])

    def wildtype_sequence(self):
        return ''.join([nucleotide[1] for nucleotides in self.sequence.values() for nucleotide in nucleotides])

    def altered_sequence(self):
        return ''.join([nucleotide[0] for nucleotides in self.sequence.values() for nucleotide in nucleotides])

    def alignment(self):
        return self.wildtype_sequence(), self.alignment_string(), self.altered_sequence()

    def seq_from_original_coordinates(self, start, end):
        """Get edited sequence from coordinates of wild type sequence"""
        start += 1
        end += 1
        seq = []
        length = end - start
        for i, pos in enumerate(list(self.sequence.values())[start:end], 1):
            for j, pair in enumerate(pos):
                if i < length or j == 0:
                    if self.alteration_list[start+i-2][j] != 'n':
                        seq.append(pair[1].lower())
                    else:
                        seq.append(pair[1])
        return ''.join(seq)

    def seq_from_new_coordinates(self, start, end):
        """Get wild type sequence from coordinates of edited sequence"""
        seq = []
        deletions = 0
        prevpos = self.index[start]
        seq.append(self[prevpos])
        for i, pos in enumerate(self.index[start + 1:end], 1):
            while deletions + prevpos[0] < pos[0] - 1:
                seq.append('-')
                deletions += 1
            deletions = 0
            if self.alteration_list[pos[0]][pos[1]] != 'n':
                seq.append(self[pos].lower())
            else:
                seq.append(self[pos])
            prevpos = pos
        return ''.join(seq)

    @property
    def number_of_alterations(self):
        return sum([0 if n == 'n' else 1 for n in self.alteration_string()])

    @property
    def number_of_insertions(self):
        return sum([1 if n == 'i' else 0 for n in self.alteration_string()])

    @property
    def number_of_deletions(self):
        return sum([1 if n == 'd' else 0 for n in self.alteration_string()])

    def add_degenerate(self, degenerate_sequence, start):
        """Add a degenerate region to the tracker

                Required arguments:
                degenerate_sequence -- degenerate sequence of the region
                start -- first nucleotide of the degenerate region
                """
        degenerate = Degeneracy(degenerate_sequence, start)
        if start > len(self.sequence):
            return
        assert all([self[i] in degenerate[i] for i in range(degenerate.start, min(len(self.index), degenerate.end))])

        self.degeneracy.append(degenerate)

    def degenerate_sequence(self, not_dgn='N', deletion_symbol='', strict=True):
        """Get a degenerate representation of the edited sequence.

                Optional arguments:
                not_dgn -- symbol to use for nucleotides not in a degenerate region, defaults to N
                deletion_symbol -- symbol to use for deleted nucleotides, defaults to no symbol
                strict -- Whether to represent nucleotides outside degenerate regions as 'free' nucleotides or 'locked'
                          nucleotides, defaults to locked.
                """
        dgn_str = []
        position_index = {idx: i for i, idx in enumerate(self.index)}

        for i, (pos, alts) in enumerate(zip(self.sequence.values(), self.alteration_list.values()), -1):
            for j, (nts, alt) in enumerate(zip(pos, alts)):
                if alt == 'd':
                    dgn_str.append(deletion_symbol)
                else:
                    position = position_index[(i, j)]
                    for dgn in self.degeneracy:
                        if position in dgn:
                            dgn_str.append(str(dgn)[position])
                            break
                    else:
                        if strict:
                            dgn_str.append(nts[0])
                        else:
                            dgn_str.append(not_dgn)
        return ''.join(dgn_str)

    @property
    def alteration_length(self):
        """Length of the edited region"""
        alterations = self.alterations()

        if not alterations:
            return 0

        insertions = self.number_of_insertions
        deletions = self.number_of_deletions

        start = sum(alterations[0])
        end = alterations[-1][0] + 1

        return max(end-start, len(self.index[start:end + insertions - deletions]))

    def make_oligos(self, repair, num_pegs, nuclease=None, silence_pam=False, **options):
        """Make oligos for all selected spacers"""
        if not self.number_of_alterations:
            raise Exception('No DNA changes found')
        print(nuclease, options)
        oligo_sets = self.find_best_spacers(num_pegs, repair, nuclease, **options)
        degenerate_sequence = self.degenerate_sequence(strict=silence_pam == 'strict')
        for oligo_set in oligo_sets:
            oligo_set.make_oligos(degenerate_sequence, silence_pam, **options)

        # Visualisations for Detail view
        alterations = self.alterations()

        # Visualisations for List view
        # TODO: FIX FOR CASES WITH ONLY SPACERS ON ONE STRAND !!!

        start = sum(alterations[0])
        end = sum(alterations[-1])

        for oligo_set in oligo_sets:
            if oligo_set.spacer_strand == 1:
                if oligo_set.spacer_position < start:
                    start = oligo_set.spacer_position

            else:
                if oligo_set.spacer_position > end:
                    end = oligo_set.spacer_position

        start -= 20
        end += 20
        start = max(0, start)
        if repair:
            end = min(end, len(self.index))
            visual_sequence = self.seq_from_new_coordinates(start, end)
            diff = 180 - len(visual_sequence)
            if diff > 0:
                start -= diff//2
                end += diff//2
                start = max(0, start)
                end = min(end, len(self.index))
            visual_sequence = self.seq_from_new_coordinates(start, end)

            nicking_start = max(start - (options['nicking_range'] + 50), 0)
            nicking_end = min(end + (options['nicking_range'] + 50), len(self.index))
            visual_nicking = self.seq_from_new_coordinates(nicking_start, nicking_end)

        else:
            end = min(end, len(self._original_sequence))
            diff = 180 - (end - start)
            if diff > 0:
                start -= diff//2
                end += diff//2
                start = max(0, start)
                end = min(end, len(self._original_sequence))
            visual_sequence = self.seq_from_original_coordinates(start, end)

            nicking_start = max(start - (options['nicking_range'] + 50), 0)
            nicking_end = min(end + (options['nicking_range'] + 50), len(self._original_sequence))
            visual_nicking = self.seq_from_original_coordinates(nicking_start, nicking_end)

        indels = self.number_of_insertions + self.number_of_deletions
        for o in oligo_sets:
            o.spacer_position -= start
            if o.spacer_strand == -1:
                o.spacer_position += visual_sequence.count('-')
            for n in o.nicking_spacers:
                n['position'] += indels

        translations = []
        upstream = ''
        offset = 0
        aa_index = 0
        for i, dgn in enumerate(self.degeneracy):
            offset += (dgn.end - dgn.start) % 3
            offset = offset % 3
            if dgn.end >= start:
                translations.append({'AA_index': aa_index, 'upstream': upstream, 'start': dgn.start-start, 'end': dgn.end-start})
            if dgn.start > end:
                break

            aa_index += (dgn.end - dgn.start) // 3
            upstream = self[dgn.end-offset:dgn.end]

        nicking_offset = visual_nicking.find(visual_sequence)

        return {
            'pegRNAs': sorted([oligo_set.info for oligo_set in oligo_sets],
                              key=lambda os: (-os['pam_disrupted'], -os['pam_silenced'], os['distance'])),
            'primers': self.make_primers(**options),
            'visual_sequence': visual_sequence,
            'visual_nicking': visual_nicking,
            'nicking_offset': nicking_offset,
            'translations': translations,
            'start': start,
        }

    def find_best_spacers(self, num_pegs, repair=False, nuclease=None, **options):
        """Find pegRNA spacers.

                Required arguments:
                num_pegs -- number of pegRNAs to return

                Optional arguments:
                repair -- If true, designs pegRNAs from edited sequence -> wild type sequence, defaults to False.
                nuclease -- Which nuclease to design pegRNAs for, if None uses default nuclease.
                """
        reference_sequence = self.original_sequence
        altered_sequence = self.__str__()

        alterations = self.alterations()
        position = sum(alterations[0])

        if nuclease is None:
            nuclease = django.conf.settings.DESIGN_CONF['default_nuclease']

        nuclease = NUCLEASES[nuclease]

        if repair:
            reference_sequence, altered_sequence = altered_sequence, reference_sequence

        spacers = nuclease.find_spacers(reference_sequence, altered_sequence, position, position + self.alteration_length, **options)
        n = min(num_pegs, len(spacers))
        for i in range(n):
            spacer = spacers[i]['spacer']
            pos = spacers[i]['position']
            strand = spacers[i]['strand']

            if strand == -1:
                pos = pos - len(spacer) + 1
            if repair:
                visual_spacer = self.seq_from_new_coordinates(pos, pos + len(spacer))
            else:
                visual_spacer = self.seq_from_original_coordinates(pos, pos + len(spacer))
            if strand == -1:
                visual_spacer = reverse_complement(visual_spacer)
            spacers[i]['visual_spacer'] = visual_spacer

        return [OligoSet(tracker=self, spacer=spacers[i], nuclease=nuclease, repair=repair)
                for i in range(n)]

    def make_primers(self, product_min_size, product_max_size, primer_min_length, primer_max_length, primer_opt_length,
                     primer_min_tm, primer_max_tm, primer_opt_tm, **options):
        """Design sequencing primers to verfiy edit."""
        alterations = self.alterations()
        start = sum(alterations[0])
        target = str(self)
        target = target[max(0, start-product_max_size):min(start+self.alteration_length - self.number_of_deletions + product_max_size, len(target))]

        if (start - product_max_size) > 0:
            start = product_max_size

        out = run_primer3(target, start, self.alteration_length - self.number_of_deletions, product_min_size,
                          product_max_size, primer_min_length, primer_max_length, primer_opt_length, primer_min_tm,
                          primer_max_tm, primer_opt_tm, **options)
        primer_attribute = regex.compile(r'PRIMER_(?P<PRIMER>\w+)_(?P<NUMBER>\d+)(?=_(?P<ATTRIBUTE>\w+))*')
        primer_sequence = regex.compile(r'PRIMER_(?P<PRIMER>\w+)_(?P<NUMBER>\d+)_SEQUENCE*')

        primers = defaultdict(lambda: {'primers': {'LEFT': {}, 'RIGHT': {}, 'PAIR': {}}, 'products': list()})

        for i, line in enumerate(out.split('\n')):
            if line.startswith('='):
                break
            if i < 17:
                continue

            label, value = line.split('=')
            m = primer_sequence.match(label)
            a = primer_attribute.match(label)
            if m:
                number = m.group('NUMBER')
                primer = m.group('PRIMER')
                primers[f'P{number}']['primers'][primer]['SEQUENCE'] = value
            elif a:
                number = a.group('NUMBER')
                primer = a.group('PRIMER')
                attribute = a.group('ATTRIBUTE')
                if attribute is None:
                    continue
                primers[f'P{number}']['primers'][primer][attribute] = value

        for pair in primers:
            primers[pair]['products'] = list(primers[pair]['products'])
        return [{**pair, 'id': id} for id, pair in primers.items()]


class Degeneracy:
    """Degenerate nucleotides for alteration tracking."""
    def __init__(self, degenerate_sequence, start):
        if not is_valid_degenerate(degenerate_sequence):
            raise ValueError(f'Sequence is not valid degenerate DNA {degenerate_sequence}')
        self.sequence = degenerate_sequence.upper()

        self.start = start
        self.end = start + len(degenerate_sequence)

        self.allowed = [degenerate_to_nucleotides[nucleotide] for nucleotide in self.sequence]

    def __getitem__(self, item):
        return self.allowed[item - self.start]

    def __str__(self):
        return ' ' * self.start + self.sequence

    def __repr__(self):
        return f'{self.start}-{self.end}'

    def __contains__(self, i):
        return self.start <= i < self.end

    def delete(self, i):
        """Delete a dgn"""
        del self.allowed[i - self.start]
        sequence = list(self.sequence)
        del sequence[i - self.start]
        self.sequence = ''.join(sequence)
        self.end -= 1

    def insert(self, i, dgn):
        """Insert a dgn"""
        self.allowed.insert(i - self.start, degenerate_to_nucleotides[dgn.upper()])
        sequence = list(self.sequence)
        sequence.insert(i - self.start, dgn.upper())
        self.sequence = ''.join(sequence)
        self.end += 1

    def substitute(self, i, dgn):
        """Make a substitution"""
        self.allowed[i - self.start] = degenerate_to_nucleotides[dgn.upper()]
        sequence = list(self.sequence)
        sequence[i - self.start] = dgn.upper()
        self.sequence = ''.join(sequence)
