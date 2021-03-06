from unittest import TestCase

from design.nucleases.Cas9 import SpCas9, Anzalone, GGAssembly, LibraryCloning


class GGAssemblyTestCase(TestCase):

    def test_make_scaffold_oligos(self):
        oligos = GGAssembly.make_scaffold_oligos('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC')
        self.assertEqual(oligos, {'top': 'AGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCG', 'bottom': 'GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAG'})

    def test_make_spacer_oligos(self):
        # Starts with G
        oligos = GGAssembly.make_spacer_oligos('GGCCCAGACTGAGCACGTGA', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC')
        self.assertEqual(oligos, {'top': 'caccGGCCCAGACTGAGCACGTGAgtttt', 'bottom': 'ctctaaaacTCACGTGCTCAGTCTGGGCC'})
        # Doesn't start with G
        oligos = GGAssembly.make_spacer_oligos('CGCCCAGACTGAGCACGTGA', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC')
        self.assertEqual(oligos, {'top': 'caccgCGCCCAGACTGAGCACGTGAgtttt', 'bottom': 'ctctaaaacTCACGTGCTCAGTCTGGGCGc'})

    def test_make_extension_oligos(self):
        oligos = GGAssembly.make_extension_oligos('TCTGCCATCAAAGCGTGCTCAGTCTG', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC')
        self.assertEqual(oligos, {'top': 'gtgcTCTGCCATCAAAGCGTGCTCAGTCTG', 'bottom': 'aaaaCAGACTGAGCACGCTTTGATGGCAGA'})

    def test_make_nicking_oligos(self):
        # Starts with G
        oligos = GGAssembly.make_nicking_oligos('GGCCCAGACTGAGCACGTGA', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC')
        self.assertEqual(oligos, {'top': 'caccGGCCCAGACTGAGCACGTGA', 'bottom': 'aaacTCACGTGCTCAGTCTGGGCC'})
        # Doesn't start with G
        oligos = GGAssembly.make_nicking_oligos('CGCCCAGACTGAGCACGTGA', 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC')
        self.assertEqual(oligos, {'top': 'caccgCGCCCAGACTGAGCACGTGA', 'bottom': 'aaacTCACGTGCTCAGTCTGGGCGc'})


class SpCas9TestCase(TestCase):
    pass


class AnzaloneTestCase(TestCase):

    def test_make_rt_sequence(self):
        reference = 'cttTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        self.assertEqual(Anzalone._make_rt_sequence(SpCas9, reference, 3, 0, 3, 10, 34, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'cttTGATGGCAGA')
        # Should be same as above, extension should not start with a C
        self.assertEqual(Anzalone._make_rt_sequence(SpCas9, reference, 3, 0, 3, 9, 34, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'cttTGATGGCAGA')
        reference = 'aTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        self.assertEqual(Anzalone._make_rt_sequence(SpCas9, reference, 1, 0, 1, 10, 34, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'aTGATGGCAGA')
        self.assertEqual(Anzalone._make_rt_sequence(SpCas9, reference, 1, 0, 1, 9, 34, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'aTGATGGCAGA')
        reference = 'GATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        self.assertEqual(Anzalone._make_rt_sequence(SpCas9, reference, 1, -1, 1, 9, 34, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'GATGGCAGA')
        self.assertEqual(Anzalone._make_rt_sequence(SpCas9, reference, 1, -1, 1, 10, 34, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'GATGGCAGAGGA')

    def test_make_pbs_sequence(self):
        reference = 'TATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACG'
        self.assertEqual(Anzalone._make_pbs_sequence(SpCas9, reference, 13, 20, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'CAGACTGAGCACG')
        self.assertEqual(Anzalone._make_pbs_sequence(SpCas9, reference, 10, 10, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'ACTGAGCACG')
        self.assertEqual(Anzalone._make_pbs_sequence(SpCas9, reference, 11, 11, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'GACTGAGCACG')
        self.assertEqual(Anzalone._make_pbs_sequence(SpCas9, reference, 10, 20, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'ACTGAGCACG')
        reference = 'TATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAAAAAA'
        self.assertEqual(Anzalone._make_pbs_sequence(SpCas9, reference, 13, 20, strategy=list(SpCas9.cloning_strategies.keys())[0])[0], 'CCCAGACTGAAAAAA')

    def test_make_extension_sequence(self):
        pbs_min_length=13
        pbs_max_length=20
        rt_min_length=10
        rt_max_length=34

        # CTT_insertion
        wt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        alt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGcttTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        pbs_length, rt_length, extension_seq, alternates = Anzalone.make_extension_sequence(SpCas9, wt, alt, 1, 200, 3, 3, pbs_min_length, pbs_max_length, rt_min_length, rt_max_length, strategy=list(SpCas9.cloning_strategies.keys())[0])
        self.assertEqual(extension_seq, 'TCTGCCATCAaagCGTGCTCAGTCTG')

        # T to A
        alt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGaGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        pbs_length, rt_length, extension_seq, alternates = Anzalone.make_extension_sequence(SpCas9, wt, alt, 1, 200, 1, 1, pbs_min_length, pbs_max_length, rt_min_length, rt_max_length, strategy=list(SpCas9.cloning_strategies.keys())[0])
        self.assertEqual(extension_seq, 'TCCTCTGCCATCtCGTGCTCAGTCTG')

        # T deletion
        alt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        pbs_length, rt_length, extension_seq, alternates = Anzalone.make_extension_sequence(SpCas9, wt, alt, 1, 200, 1, 1, pbs_min_length, pbs_max_length, rt_min_length, rt_max_length, strategy=list(SpCas9.cloning_strategies.keys())[0])
        self.assertEqual(extension_seq, 'TCCTCTGCCATCCGTGCTCAGTCTG')

        # A insertion
        alt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGaTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        pbs_length, rt_length, extension_seq, alternates = Anzalone.make_extension_sequence(SpCas9, wt, alt, 1, 200, 1, 1, pbs_min_length, pbs_max_length, rt_min_length, rt_max_length, strategy=list(SpCas9.cloning_strategies.keys())[0])
        self.assertEqual(extension_seq, 'TCTGCCATCAtCGTGCTCAGTCTG')

    def test_find_spacers(self):
        wt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        alt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGcttTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        spacers = Anzalone.find_spacers(SpCas9, wt, alt, 200, 203, 100)
        self.assertEqual(spacers[0]['spacer'], 'GGCCCAGACTGAGCACGTGA')
        self.assertEqual(spacers[0]['position'], 183)
        self.assertEqual(spacers[0]['cut_site'], 200)
        self.assertEqual(spacers[0]['strand'], 1)
        self.assertEqual(spacers[0]['pam'], ('TGG', 203))
        self.assertEqual(spacers[0]['pam_disrupted'], True)
        self.assertEqual(spacers[0]['distance'], 0)
        alt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGaGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        spacers = Anzalone.find_spacers(SpCas9, wt, alt, 200, 201, 100)
        self.assertEqual(spacers[0]['spacer'], 'GGCCCAGACTGAGCACGTGA')
        self.assertEqual(spacers[0]['position'], 183)
        self.assertEqual(spacers[0]['cut_site'], 200)
        self.assertEqual(spacers[0]['strand'], 1)
        self.assertEqual(spacers[0]['pam'], ('TGG', 203))
        self.assertEqual(spacers[0]['pam_disrupted'], False)
        self.assertEqual(spacers[0]['distance'], 0)

    def test_find_nicking_spacers(self):
        wt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        alt = 'TGGAGAGTTTTAAGCAAGGGCTGATGTGGGCTGCCTAGAAAGGCATGGATGAGAGAAGCCTGGAGACAGGGATCCCAGGGAAACGCCCATGCAATTAGTCTATTTCTGCTGCAAGTAAGCATGCATTTGTAGGCTTGATGCTTTTTTTCTGCTTCTCCAGCCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGcttTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGCGTCGCAGGACAGCTTTTCCTAGACAGGGGCTAGTATGTGCAGCTCCTGCACCGGGATACTGGTTGACAAGTTTGGCTGGGCTGGAAGCCAGCACCTAGGGAGGTCCCTGGAAGGGGCCAGCCTCACCAGGAGAGGAGGGACCTGGCCCTTCAGGGTCGAG'
        spacers = Anzalone.find_nicking_spacers(SpCas9, wt, alt, 1, 200, 'sgRNA', 100, GGAssembly)

        self.assertEqual(spacers[0]['spacer'], 'CCATCAAAGCGTGCTCAGTC')
        self.assertEqual(spacers[0]['position'], -8)
        self.assertEqual(spacers[0]['kind'], '3b')

        self.assertEqual(spacers[2]['kind'], '3')
        self.assertEqual(spacers[2]['spacer'], 'GCACATACTAGCCCCTGTCT')
        self.assertEqual(spacers[2]['position'], 66)