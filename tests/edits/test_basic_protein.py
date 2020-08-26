from unittest import TestCase

from design.edits.basic_protein import AminoAcidAlteration, TagProtein
from design.helpers import DEFAULT_CODON_USAGE_TABLE, translate


class SequenceObject:
    """Test sequence object."""
    def __init__(self, sequence, coding_sequence=None):
        self.sequence = sequence
        self.upstream = ''
        self.downstream = ''
        self.coding_regions = coding_sequence
        self.codon_usage_table = DEFAULT_CODON_USAGE_TABLE
        if coding_sequence is None:
            self.coding_sequence = ''
        else:
            cds = []
            for (start, end) in coding_sequence:
                cds.append(sequence[start:end])
            self.coding_sequence = ''.join(cds)


    @staticmethod
    def contextualize_alteration(alteration):
        return alteration

    @property
    def translation(self):
        return translate(self.coding_sequence)

    @property
    def degenerate_sequence(self):
        return []

    def cds_pos_to_position(self, position):
        distance = 0
        for region in self.coding_regions:
            start, end = region

            length = end - start + 1

            if position < distance + length:
                return start + position - distance
            distance += length
        return None

    def protein_to_transcript_position(self, protein_position, insertions_before=0):
        protein_position = protein_position - insertions_before
        start, end = protein_position * 3, protein_position * 3 + 2

        positions = [self.cds_pos_to_position(position) for position in list(range(start, end + 1))]
        return [p+insertions_before*3 for p in positions]


class BasicProteinTestCase(TestCase):

    def setUp(self) -> None:
        self.sequence_object = SequenceObject('TTTTTTTTTTATGCGAATATTGTTTTTTTTTT', [(10, 22)])

    def test_amino_acid_alteration(self):
        edit = AminoAcidAlteration(self.sequence_object, option_string='alteration=R2A')
        tracker = edit.run()
        self.assertEqual(str(tracker), 'TTTTTTTTTTATGgcAATATTGTTTTTTTTTT')

    def test_tag_protein(self):
        edit = TagProtein(self.sequence_object, option_string='tag=FLAG,position=N')
        tracker = edit.run()
        self.assertEqual(str(tracker), 'TTTTTTTTTTATGgactacaaggacgacgacgacaagCGAATATTGTTTTTTTTTT')
