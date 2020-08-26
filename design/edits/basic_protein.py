"""
Module for basic protein alterations
"""
import itertools
import regex

from . import AbstractEdit, register_edit

from design.helpers import DEFAULT_CODON_TABLE, DEGENERATE_CODON_TABLE, translate, degenerate_to_nucleotides

amino_acid_alteration_format = regex.compile(r'^([ACDEFGHIKLMNPQRSTVWY*$]+)(\d+)([ACDEFGHIKLMNPQRSTVWY*$]+)$')


class AminoAcidAlteration(AbstractEdit):
    """Substitute amino acids.

    ###### Options:
    * *alteration (str)*: List of amino acid substitions.
    Multiple substitutions can be separated by a semicolon(;).
    * Deletions are denoted with an asterisk (\*)
        * 5N\*.
    * Insertions are denoted with an asterisk (\*)
        * \*5A.
    * Stop codons are denoted with a dollar sign ($)
        * N5$.

    ###### Usage
    Protein sequence is: MVHLTPEEK
    * alteration=H3A
        * New protein sequence is: MVALTPEEK

    * alteration=H3*
        * New protein sequence is: MVLTPEEK

    * alteration=*3A
        * New protein sequence is: MVAHLTPEEK

    Neighbouring alterations can be changed in one substitution
    * alteration=HL3AR
        * New protein sequence is: MVARTPEEK

    * alteration=HL3*
        * New protein sequence is: MVTPEEK

    * alteration=*3AR
        * New protein sequence is: MVARHLTPEEK

    """

    default_options = {'alteration': (str, None, "^(([ACDEFGHIKLMNPQRSTVWY*$]+)(\d+)([ACDEFGHIKLMNPQRSTVWY*$]+);)*(([ACDEFGHIKLMNPQRSTVWY*$]+)(\d+)([ACDEFGHIKLMNPQRSTVWY*$]+);?)$")}
    required = ('alteration',)
    allow_genomic = False

    @classmethod
    def validate_options(cls, options, sequence_object):
        super().validate_options(options, sequence_object)
        alterations = options['alteration'].split(';')
        parsed_alterations = []
        codon_table = DEFAULT_CODON_TABLE

        for alteration in alterations:
            """"""
            if alteration:
                m = regex.match(amino_acid_alteration_format, alteration)
                if not m:
                    raise ValueError(f"Invalid alteration '{alteration}' of format 'A10B' in mut={options['alteration']}")

                old_amino_acid, position, new_amino_acid = m.groups()
                position = int(position)

                for i, (o, n) in enumerate(itertools.zip_longest(old_amino_acid, new_amino_acid, fillvalue='*')):
                    parsed_alterations.append((o, position + i, n))
                protein_length = int(len(sequence_object.coding_sequence) / 3)
                for (old_amino_acid, position, new_amino_acid) in parsed_alterations:

                    gene_position = (position - 1) * 3
                    if old_amino_acid == new_amino_acid:
                        raise ValueError(f"Invalid amino acid change: {''.join([old_amino_acid, str(position), new_amino_acid])}.")
                    if position > protein_length:
                        raise ValueError(f"Invalid amino acid change: {''.join([old_amino_acid, str(position), new_amino_acid])}. Protein sequence is {protein_length} AA long.")
                    expected_codon = sequence_object.coding_sequence[gene_position:gene_position + 3]
                    if old_amino_acid != '*':
                        expected_amino_acid = codon_table[expected_codon]
                        if str(old_amino_acid) != str(expected_amino_acid):
                            raise ValueError(f"Invalid amino acid change: {''.join([old_amino_acid, str(position), new_amino_acid])}. Expecting to change {old_amino_acid}, but {expected_amino_acid} was found in the translation.")
        return parsed_alterations, codon_table

    def run(self):
        """"""
        parsed_alterations, codon_table = self.validate_options(self.options, self.sequence_object)
        return self.run_parsed(parsed_alterations, codon_table)

    def run_parsed(self, parsed_alterations, codon_table):
        degenerate_codon_table = DEGENERATE_CODON_TABLE
        codon_usage_table = self.sequence_object.codon_usage_table

        offset = 0
        insertions = 0
        tracker = self.get_tracker()
        max_dna_pos = max(self.sequence_object.protein_to_transcript_position(len(self.sequence_object.translation)-1))
        for (old_amino_acid, position, new_amino_acid) in parsed_alterations:

            gene_position = (position - 1) * 3
            expected_codon = self.sequence_object.coding_sequence[gene_position:gene_position + 3]
            if old_amino_acid != '*':
                expected_amino_acid = codon_table[expected_codon]
                if str(old_amino_acid) != str(expected_amino_acid):
                    raise ValueError(f"Invalid amino acid change: {''.join([old_amino_acid, str(position), new_amino_acid])}. Expecting to change {old_amino_acid}, but {expected_amino_acid} was found in the translation.")

            codon_positions = self.sequence_object.protein_to_transcript_position(position - 1 + offset, insertions)
            for i, codon_pos in enumerate(codon_positions):
                if codon_pos is None:
                    max_dna_pos += 1
                    codon_positions[i] = max_dna_pos
            codon_positions = [self.sequence_object.contextualize_alteration(position) for position in codon_positions]
            possible_codons = [codon for codon, residue in codon_table.items() if residue == new_amino_acid]

            if new_amino_acid == '*':
                positions = list(itertools.chain.from_iterable(
                    [self.sequence_object.protein_to_transcript_position(position + i + offset) for i in range(-2, 1)]))
                positions_to_delete = self.amino_acid_positions_to_delete(positions, str(tracker),
                                                                          degenerate_codon_table)
                for i, pos in enumerate(positions_to_delete):
                    pos = self.sequence_object.contextualize_alteration(pos)
                    tracker.delete(pos-i, force=True)
                offset -= 1

            elif old_amino_acid == '*':
                codon = max([(codon, codon_usage_table[codon][1]) for codon in possible_codons], key=lambda x: x[1])[0]
                for i, nucleotide in enumerate(codon):
                    tracker.insert(nucleotide, codon_positions[0]+i, force=True)
                insertions += 1
            else:
                new_codon = possible_codons[0]
                changed_positions = [codon_positions[i] for i, (n, o) in enumerate(zip(new_codon, expected_codon)) if
                                     n != o]
                minium_distance = max(changed_positions) - min(changed_positions)
                for codon in possible_codons[1:]:
                    changed_positions = [codon_positions[i] for i, (n, o) in enumerate(zip(codon, expected_codon)) if n != o]
                    distance = max(changed_positions) - min(changed_positions)
                    # Distance is longer than best, skip
                    if distance > minium_distance:
                        continue
                    # Distance is shorter than best, keep
                    if distance < minium_distance:
                        minium_distance = distance
                        new_codon = codon
                    # Distance is the same, keep if usage is better
                    else:
                        if codon_usage_table[codon] > codon_usage_table[new_codon]:
                            new_codon = codon

                for i, (pos, nucleotide) in enumerate(zip(codon_positions, new_codon)):
                    tracker.mutate(nucleotide, pos, force=True)

        return tracker


    @staticmethod
    def amino_acid_positions_to_delete(positions, sequence, degenerate_codon_table):

        degenerate_upstream = degenerate_codon_table[translate(''.join(sequence[pos] for pos in positions[:3]))]
        degenerate_downstream = degenerate_codon_table[translate(''.join(sequence[pos] for pos in positions[6:]))]
        codon = ''.join(sequence[pos] for pos in positions[3:6])

        # If all nucleotides in codon to delete are on same exon, delete codon
        if positions[5] - positions[3] == 2:
            return positions[3:6]
        # First two nucleotides in same exon
        if positions[4] - positions[3] == 1:
            # Can last nucleotide of codon be used as last nucleotide of upstream codon?
            if degenerate_to_nucleotides[codon[2]] <= degenerate_to_nucleotides[degenerate_upstream[2]]:
                return positions[2:5]
            # Can first two nucleotides of codon be used as first two nucleotides of downstream codon?
            if all([degenerate_to_nucleotides[c] <= degenerate_to_nucleotides[d] for (c, d) in zip(codon[:2],
                                                                                                   degenerate_downstream[
                                                                                                   :2])]):  # if degenerate_to_nucleotides[codon[0]] <= degenerate_to_nucleotides[degenerate_downstream[0]] and degenerate_to_nucleotides[codon[1]] <= degenerate_to_nucleotides[degenerate_downstream[1]]:
                return positions[5:8]

        # Last two nucleotides in same exon
        elif positions[5] - positions[4] == 1:
            # Can last two nucleotides of codon be used as last two nucleotides of upstream codon?
            if all([degenerate_to_nucleotides[c] <= degenerate_to_nucleotides[u] for (c, u) in zip(codon[1:],
                                                                                                   degenerate_upstream[
                                                                                                   1:])]):  # if degenerate_to_nucleotides[codon[1]] <= degenerate_to_nucleotides[degenerate_upstream[1]] and degenerate_to_nucleotides[codon[2]] <= degenerate_to_nucleotides[degenerate_upstream[2]]:
                return positions[1:4]
            # Can first nucleotide of codon be used as first nucleotide of downstream codon?
            if degenerate_to_nucleotides[codon[0]] <= degenerate_to_nucleotides[degenerate_downstream[0]]:
                return positions[4:7]

        # Must delete nucleotides in both exons
        return positions[3:6]


class TagProtein(AminoAcidAlteration):
    """Add a tag to a protein.

    ###### Options:
    * *tag (select)*: Select a tag from the list.
    * *position (int / C / N)*: Position of the tag. N-terminal or C-terminal or position of insertion.

    Available tags:
    * 5xHIS: HHHHH,
    * 6xHIS: HHHHHH,
    * FLAG: DYKDDDDK,
    * HA: YPYDVPDYA,
    * HSV: QPELAPEDPED,
    * Myc: EQKLISEEDL,
    * V5: GKPIPNPLLGLDST,
    """

    tags = {
        '5xHIS': 'HHHHH',
        '6xHIS': 'HHHHHH',
        'FLAG': 'DYKDDDDK',
        'HA': 'YPYDVPDYA',
        'HSV': 'QPELAPEDPED',
        'Myc': 'EQKLISEEDL',
        'V5': 'GKPIPNPLLGLDST',
    }

    default_options = {'tag': (str, None, list(tags.keys())),
                       'position': (str, None, '[CcNn]+|\d+')}
    required = ('tag', 'position')

    @classmethod
    def validate_options(cls, options, sequence_object):
        super(AminoAcidAlteration, cls).validate_options(options, sequence_object)
        tag = options['tag']

        position = options['position']
        try:
            position = int(position)
        except ValueError:
            if position.upper() == 'N':
                position = 1
                if sequence_object.translation[0] == 'M':
                    position = 2
            elif position.upper() == 'C':
                position = len(sequence_object.translation)
            else:
                raise ValueError(f"Invalid position {position}. Allowed values are, C, N or a number.")

        if tag not in cls.tags.keys():
            raise ValueError(f"Invalid tag: {tag}. Tag mut be one of: '{','.join(cls.tags.keys())}'")

        return tag, position

    def run(self):
        """"""
        tag, position = self.validate_options(self.options, self.sequence_object)
        codon_table = DEFAULT_CODON_TABLE
        parsed_alterations = []
        for i, amino_acid in enumerate(self.tags[tag]):
            parsed_alterations.append(('*', position + i, amino_acid))

        return self.run_parsed(parsed_alterations, codon_table)


register_edit(AminoAcidAlteration)
register_edit(TagProtein)
