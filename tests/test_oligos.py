from unittest import TestCase

from design.oligos import AlterationTracker


class AlterationTrackerTestCase(TestCase):

    def setUp(self) -> None:
        self.tracker = AlterationTracker('ATGCAGCTAGCAGTACA')

    def test_insert(self):
        self.tracker.insert('A', 2)
        self.assertEqual(str(self.tracker), 'ATaGCAGCTAGCAGTACA')
        self.tracker.insert('T', 2)
        self.assertEqual(str(self.tracker), 'ATtaGCAGCTAGCAGTACA')
        self.assertEqual(self.tracker.alteration_string(), 'nniinnnnnnnnnnnnnnn')
        self.assertEqual(self.tracker.number_of_insertions, 2)
        self.assertEqual(self.tracker.number_of_deletions, 0)
        self.assertEqual(self.tracker.number_of_alterations, 2)

    def test_delete(self):
        self.tracker.delete(2)
        self.assertEqual(str(self.tracker), 'ATCAGCTAGCAGTACA')
        self.tracker.delete(3)
        self.assertEqual(str(self.tracker), 'ATCGCTAGCAGTACA')
        self.assertEqual(self.tracker.alteration_string(), 'nndndnnnnnnnnnnnn')
        self.assertEqual(self.tracker.number_of_insertions, 0)
        self.assertEqual(self.tracker.number_of_deletions, 2)
        self.assertEqual(self.tracker.number_of_alterations, 2)

    def test_substitution(self):
        self.tracker.substitute('A', 2)
        self.assertEqual(str(self.tracker), 'ATaCAGCTAGCAGTACA')
        self.tracker.substitute('A', 5)
        self.assertEqual(str(self.tracker), 'ATaCAaCTAGCAGTACA')
        self.assertEqual(self.tracker.alteration_string(), 'nnmnnmnnnnnnnnnnn')
        self.assertEqual(self.tracker.number_of_insertions, 0)
        self.assertEqual(self.tracker.number_of_deletions, 0)
        self.assertEqual(self.tracker.number_of_alterations, 2)
