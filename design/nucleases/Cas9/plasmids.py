import random

from design.nucleases import BaseCloningStrategy
from design.nucleases import OligoDict
from design.helpers import reverse_complement


class GGAssembly(BaseCloningStrategy):
    """Cloning as in Anzalone et al, 2019."""

    excel_oligo_headers = ['spacer_oligo_top', 'spacer_oligo_bottom', 'extension_oligo_top', 'extension_oligo_bottom']
    excel_extension_headers = ['top_oligo', 'bottom_oligo']

    @classmethod
    def help_text(cls, scaffold):
        scaffold = cls.make_scaffold_oligos(scaffold)
        return f"""Cloning as in Anzalone et al, 2019.

        Golden Gate Assembly of spacer oligos, scaffold oligos, 
        and extension oligos into [pegRNA-GG-acceptor vector (Addgene #132777)](https://www.addgene.org/132777/).

        Scaffold oligos:
        * Top: {scaffold['top']}
        * Bottom: {scaffold['bottom']}"""

    @classmethod
    def _spacer_to_cloning(cls, spacer_sequence: str) -> str:
        spacer_sequence = spacer_sequence.upper()
        if not spacer_sequence.startswith('G'):
            spacer_sequence = 'g' + spacer_sequence

        return spacer_sequence

    @classmethod
    def make_scaffold_oligos(cls, scaffold: str) -> OligoDict:
        return {
            'top': scaffold[5:-4],
            'bottom': reverse_complement(scaffold[9:])
        }

    @classmethod
    def make_spacer_oligos(cls, spacer_sequence: str, scaffold: str) -> OligoDict:
        target = cls._spacer_to_cloning(spacer_sequence)
        return {
            'top': ''.join(['cacc', target, scaffold[:5].lower()]),
            'bottom': reverse_complement(target + scaffold[:9].lower())
        }

    @classmethod
    def make_extension_oligos(cls, extension_sequence: str, scaffold: str) -> OligoDict:
        return {
            'top': scaffold[-4:].lower() + extension_sequence,
            'bottom': reverse_complement(''.join([extension_sequence, 'tttt']))
        }

    @classmethod
    def design_cloning(cls, spacer_sequence: str, scaffold: str, extension_sequence: str, **options):
        return {
            'spacer': cls.make_spacer_oligos(spacer_sequence, scaffold),
            'scaffold': cls.make_scaffold_oligos(scaffold),
            'extension': cls.make_extension_oligos(extension_sequence, scaffold),
        }

    @classmethod
    def make_nicking_oligos(cls, spacer_sequence: str, scaffold: str):
        target = cls._spacer_to_cloning(spacer_sequence)
        return {
            'top': 'cacc' + target,
            'bottom': reverse_complement(target + scaffold[:4].lower())
        }

    @classmethod
    def alternate_extension(cls, extension_sequence: str, scaffold: str, **options):
        return cls.make_extension_oligos(extension_sequence, scaffold)


class LibraryCloning(BaseCloningStrategy):
    """Cloning for libraries."""

    can_design_primers = False
    can_design_nicking = False
    excel_oligo_headers = ['assembly_oligo']
    excel_extension_headers = ['assembly_oligo']
    allow_extension_filtering = False

    @classmethod
    def help_text(cls, scaffold):
        scaffold = cls.make_scaffold_oligos(scaffold)
        return f"""Cloning into modified lentiGuide-Puro plasmid for construction of prime editing libraries.

    Gibson/NeBuilder assembly into modified lentiGuide-puro vector.

    5’- atcttgtggaaaggacgaaacacc - pegRNA_spacer - SCAFFOLD CLONING SITE – pegRNA extension – TTTTTTT – 15 nt barcode – TARGET SEQUENCE – aagcttggcgtaactagatc - 3’.

    The Scaffold cloning site consists of a unique 20mer barcode flanked by BsmBI restriction sites for inserting the sgRNA scaffold.

    No primers or nicking sgRNAs are designed when using this cloning method.

    Oligos for cloning the scaffold into assembled plasmids:
    * Top: {scaffold['top']}
    * Bottom: {scaffold['bottom']}"""

    @classmethod
    def make_scaffold_oligos(cls, scaffold: str) -> OligoDict:
        return {
            'top': scaffold[5:-4],
            'bottom': reverse_complement(scaffold[9:])
        }

    @classmethod
    def _spacer_to_cloning(cls, spacer_sequence: str) -> str:
        spacer_sequence = spacer_sequence.upper()
        if not spacer_sequence.startswith('G'):
            spacer_sequence = 'g' + spacer_sequence

        return spacer_sequence

    @classmethod
    def design_cloning(cls, spacer_sequence: str, scaffold: str, extension_sequence: str, **options):
        upstream = 'atcttgtggaaaggacgaaacacc'
        downstream = 'aagcttggcgtaactagatc'
        barcode20 = 'YYYYYYYYYYYYYYYYYYYY'
        barcode15 = 'XXXXXXXXXXXXXXX'
        target = ''.join([options['upstream'][-20:], options['downstream'][:options['cut_dist'] + 30]])
        return {'assembly': ''.join([
            upstream,
            spacer_sequence,
            scaffold[:10].lower(),
            'gagacg',  # restriction site
            barcode20,
            'cgtctc',  # restriction site
            scaffold[-5:].lower(),
            extension_sequence,
            'TTTTTTT',
            barcode15,
            target,
            downstream
        ])}

    @classmethod
    def alternate_extension(cls, spacer_sequence: str, scaffold: str, extension_sequence: str, **options):
        return cls.design_cloning(spacer_sequence, scaffold, extension_sequence, **options)

    @classmethod
    def generate_n_mer(cls, n):
        nmer = 'gagacg cgtctc'
        while any(e in nmer for e in ['gagacg', 'cgtctc']):
            nmer = ''.join([random.choice('acgt') for _ in range(n)])
        return nmer

    @classmethod
    def _update_barcodes(cls, oligo, used_15mers, used_20mers):
        bc15 = cls.generate_n_mer(15)
        bc20 = cls.generate_n_mer(20)
        while bc15 in used_15mers:
            bc15 = cls.generate_n_mer(15)
        while bc20 in used_20mers:
            bc20 = cls.generate_n_mer(20)

        oligo = oligo.replace('XXXXXXXXXXXXXXX', bc15).replace('YYYYYYYYYYYYYYYYYYYY', bc20)
        used_15mers.add(bc15)
        used_20mers.add(bc20)
        return oligo

    @classmethod
    def post_process(cls, designs):
        used_15mers = set()
        used_20mers = set()

        for design in designs:
            pegRNAs = design['pegRNAs']
            for pegRNA in pegRNAs:
                primary_oligo = pegRNA['oligos']['assembly']
                pegRNA['oligos']['assembly'] = cls._update_barcodes(primary_oligo, used_15mers, used_20mers)
                for alt_ext in pegRNA['alternate_extensions']:
                    oligo = alt_ext['oligos']['assembly']
                    alt_ext['oligos']['assembly'] = cls._update_barcodes(oligo, used_15mers, used_20mers)

            yield design
