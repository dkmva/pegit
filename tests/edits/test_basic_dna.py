from unittest import TestCase

from design.edits.basic_dna import Insertion, Deletion, Substitution


class SequenceObject:
    """Test sequence object."""
    def __init__(self, sequence):
        self.sequence = sequence
        self.upstream = ''
        self.downstream = ''

    @staticmethod
    def contextualize_alteration(alteration):
        return alteration

    @property
    def degenerate_sequence(self):
        return []


class BasicDNATestCase(TestCase):

    def setUp(self) -> None:
        self.sequence_object = SequenceObject('ATGCGAATATTG')

    def test_insertion(self):
        edit = Insertion(self.sequence_object, option_string='insert=ATA,position=7')
        tracker = edit.run()
        self.assertEqual(str(tracker), 'ATGCGAataATATTG')

    def test_deletion(self):
        edit = Deletion(self.sequence_object, option_string='delete=ATA,position=7')
        tracker = edit.run()
        self.assertEqual(str(tracker), 'ATGCGATTG')

        edit = Deletion(self.sequence_object, option_string='delete=3,position=7')
        tracker = edit.run()
        self.assertEqual(str(tracker), 'ATGCGATTG')

    def test_substitution(self):
        edit = Substitution(self.sequence_object, option_string='from=ATA,to=GGC,position=7')
        tracker = edit.run()
        self.assertEqual(str(tracker), 'ATGCGAggcTTG')