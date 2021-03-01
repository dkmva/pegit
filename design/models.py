import typing

import django 
from django.contrib.contenttypes.models import ContentType
from django.db import models, transaction
from pathlib import Path

from .helpers import parse_gff_file, extract_sequence, translate, DEGENERATE_CODON_TABLE, load_codon_table


SEQUENCE_PADDING = 600


class Organism(models.Model):
    """Organisms within design are represented by this model.

    name -- the name of an organism, eg. Homo Sapiens.
    assembly -- the name of the 2bit file containing the genomic sequence
    annotation -- the name of the annotation used for this organism, eg. a RefSeq or Gencode Release
    source -- the url where the annotation can be found
    codon_table -- the filename of a json file containing the codon usage table.
        If None, the DEFAULT_CODON_USAGE_TABLE from design.helpers will be used.
        A json codon usage table file can be writting using dump_codon_table from design.helpers
    sequence_search -- the base url for searching for a gene or transcript.
    scaffolds - text representation of a list of scaffolds in the assembly
    """
    name = models.CharField(max_length=200)
    assembly = models.CharField(max_length=200)
    annotation = models.CharField(max_length=200)
    source = models.CharField(max_length=200)
    codon_table = models.CharField(max_length=200, null=True)
    sequence_search = models.CharField(max_length=200)
    scaffolds = models.TextField(default='')

    class Meta:
        app_label = 'design'
        ordering = ('id', )

    def __str__(self):
        return f'{self.name} - {self.annotation}'

    @classmethod
    def add_to_database(cls, name, assembly, file_path, source, kind, sequence_search, codon_table=None, annotation_name=None, scaffolds='', silent=False):
        """Add an organism to the database.

        Arguments:
        name -- the name of the organism
        assembly -- name of the .2bit file containing the DNA sequences of the organism
        file_path -- path to the gff3 file containing the annotations used to populate organism gene and transcripts
        source -- source url
        kind -- ensembl or refseq - how to import the gff3 file, currently only these formats are supported.
        codon_table -- base name of the codon usage table file
        sequence_search -- sequence_search url
        annotation_name -- name of the annotation, if omitted, the file_path will be used to infer the annotation_name.
        scaffolds - text representation of a list of scaffolds in the assembly
        """
        # Dirty substitute for logging ..
        def myprint(toprint):
            if not silent:
                print(toprint)
            else:
                return
            
        batch_size = 10000
        if django.conf.settings.DATABASES['default']['ENGINE'] == 'django.db.backends.sqlite3':
            batch_size = 500

        if kind not in ['gencode', 'refseq', 'ensembl']:
            raise ValueError('Kind must be either gencode, ensembl or refseq')
        if annotation_name is None:
            annotation_name = Path(file_path).stem
        with transaction.atomic():
            organism = cls.objects.create(name=name, assembly=assembly, source=source, annotation=annotation_name,
                                          codon_table=codon_table, sequence_search=sequence_search, scaffolds=scaffolds)
            genes = []
            transcripts = []
            exons = []
            coding_sequences = []
            if kind == 'gencode':
                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']

                    if record['type'] == 'gene':
                        genes.append(Gene(
                            organism=organism,
                            gene_id=attributes['gene_id'],
                            name=attributes['gene_name'],
                            chromosome=record['chromosome'],
                            strand=record['strand'],
                            start=record['start'],
                            end=record['end'],
                            gene_type=attributes['gene_type'],
                            source=record['source']
                        ))
                myprint('insert genes')
                genes = Gene.objects.bulk_create(genes, batch_size=batch_size)
                genes = Gene.objects.filter(gene_id__in=[g.gene_id for g in genes], organism=organism)
                genes = {gene.gene_id: gene for gene in genes}
                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']
                    if record['type'] == 'transcript':
                        transcripts.append(Transcript(
                            gene=genes[attributes['gene_id']],
                            transcript_id=attributes['transcript_id'],
                            name=attributes['transcript_name'],
                            start=record['start'],
                            end=record['end'],
                            transcript_type=attributes['transcript_type'],
                            source=record['source']
                        ))
                myprint('insert transcripts')
                transcripts = Transcript.objects.bulk_create(transcripts, batch_size=batch_size)
                transcripts = Transcript.objects.filter(transcript_id__in=[t.transcript_id for t in transcripts],
                                                        gene__organism=organism)
                transcripts = {transcript.transcript_id: transcript for transcript in transcripts}
                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']
                    if record['type'] == 'exon':
                        exons.append(Exon(
                            transcript=transcripts[attributes['transcript_id']],
                            start=record['start'],
                            end=record['end'],
                        ))
                    if record['type'] == 'CDS':
                        coding_sequences.append(CodingSequence(
                            transcript=transcripts[attributes['transcript_id']],
                            start=record['start'],
                            end=record['end'],
                        ))

                myprint('insert exons')
                Exon.objects.bulk_create(exons, batch_size=batch_size)
                myprint('insert coding_sequences')
                CodingSequence.objects.bulk_create(coding_sequences, batch_size=batch_size)

            elif kind == 'refseq':
                miRNA_to_exon = {}
                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']

                    if record['type'] in ['miRNA', 'snoRNA']:
                        miRNA_to_exon[attributes['ID']] = attributes['Parent']

                    if record['type'] in ['gene', 'pseudogene']:
                        genes.append(Gene(
                            organism=organism,
                            gene_id=attributes['ID'],
                            name=attributes['Name'],
                            chromosome=record['chromosome'],
                            strand=record['strand'],
                            start=record['start'],
                            end=record['end'],
                            gene_type=attributes['gene_biotype'],
                            source=record['source']
                        ))

                myprint('insert genes')
                genes = Gene.objects.bulk_create(genes, batch_size=batch_size)
                genes = Gene.objects.filter(gene_id__in=[g.gene_id for g in genes], organism=organism)
                genes = {gene.gene_id: gene for gene in genes}

                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']
                    if record['type'] in ['transcript', 'primary_transcript', 'mRNA', 'lnc_RNA', 'antisense_RNA', 'snRNA',
                                          'guide_RNA', 'tRNA', 'V_gene_segment', 'rRNA', 'C_gene_segment', 'D_gene_segment',
                                          'V_gene_segment', 'J_gene_segment', 'telomerase_RNA', 'vault_RNA', 'Y_RNA', 'RNase_MRP_RNA',
                                          'scRNA', 'RNase_P_RNA']:
                        try:
                            name = attributes['Name']
                        except KeyError:
                            name = attributes['ID']
                        transcripts.append(Transcript(
                            gene=genes[attributes['Parent']],
                            transcript_id=attributes['ID'],
                            name=name,
                            start=record['start'],
                            end=record['end'],
                            transcript_type=genes[attributes['Parent']].gene_type,
                            source=record['source']
                        ))

                myprint('insert transcripts')
                transcripts = Transcript.objects.bulk_create(transcripts, batch_size=batch_size)
                transcripts = Transcript.objects.filter(transcript_id__in=[t.transcript_id for t in transcripts],
                                                        gene__organism=organism)
                transcripts = {transcript.transcript_id: transcript for transcript in transcripts}

                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']
                    if record['type'] == 'exon':
                        parent = attributes['Parent']
                        if parent in miRNA_to_exon.keys():
                            parent = miRNA_to_exon[parent]
                        if parent.startswith('gene-'):
                            continue
                        exons.append(Exon(
                            transcript=transcripts[parent],
                            start=record['start'],
                            end=record['end'],
                        ))
                    if record['type'] == 'CDS':
                        parent = attributes['Parent']
                        if parent in miRNA_to_exon.keys():
                            parent = miRNA_to_exon[parent]
                        if parent in genes:
                            transcripts[parent] = Transcript.objects.create(
                                gene=genes[attributes['Parent']],
                                transcript_id=attributes['Parent'],
                                name=name,
                                start=record['start'],
                                end=record['end'],
                                transcript_type=genes[parent].gene_type,
                                source=record['source']
                            )
                        coding_sequences.append(CodingSequence(
                            transcript=transcripts[parent],
                            start=record['start'],
                            end=record['end'],
                        ))
                myprint('insert exons')
                Exon.objects.bulk_create(exons, batch_size=batch_size)
                myprint('insert coding_sequences')
                CodingSequence.objects.bulk_create(coding_sequences, batch_size=batch_size)

            elif kind == 'ensembl':
                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']

                    #if record['type'] in ['ncRNA_gene', 'gene', 'pseudogene', 'transposable_element_gene']:
                    if 'ID' in attributes and attributes['ID'].split(':')[0] == 'gene':
                        try:
                            name = attributes['Name']
                        except KeyError:
                            name = attributes['ID'].split(':')[1]
                        finally:
                            genes.append(Gene(
                                organism=organism,
                                gene_id=attributes['ID'].split(':')[1],
                                name=name,
                                chromosome=record['chromosome'],
                                strand=record['strand'],
                                start=record['start'],
                                end=record['end'],
                                gene_type=attributes['biotype'],
                                source=record['source']
                            ))

                myprint('insert genes')
                genes = Gene.objects.bulk_create(genes, batch_size=batch_size)
                genes = Gene.objects.filter(gene_id__in=[g.gene_id for g in genes], organism=organism)
                genes = {gene.gene_id: gene for gene in genes}

                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']

                    #if record['type'] in ['pseudogenic_transcript', 'scRNA', 'miRNA', 'snRNA', 'tRNA', 'ncRNA', 'rRNA',
                    #                      'snoRNA', 'lnc_RNA', 'mRNA', 'transposable_element', 'pre_miRNA', 'SRP_RNA',
                    #                      'RNase_MRP_RNA']:
                    if 'ID' in attributes and attributes['ID'].split(':')[0] == 'transcript':
                        try:
                            name = attributes['Name']
                        except KeyError:
                            name = attributes['ID'].split(':')[1]
                        transcripts.append(Transcript(
                            gene=genes[attributes['Parent'].split(':')[1]],
                            transcript_id=attributes['ID'].split(':')[1],
                            name=name,
                            start=record['start'],
                            end=record['end'],
                            transcript_type=genes[attributes['Parent'].split(':')[1]].gene_type,
                            source=record['source']
                        ))

                myprint('insert transcripts')
                transcripts = Transcript.objects.bulk_create(transcripts, batch_size=batch_size)
                transcripts = Transcript.objects.filter(transcript_id__in=[t.transcript_id for t in transcripts],
                                                        gene__organism=organism)
                transcripts = {transcript.transcript_id: transcript for transcript in transcripts}

                for record in parse_gff_file(file_path, disable=silent):
                    attributes = record['attributes']
                    if record['type'] == 'exon':
                        parent = attributes['Parent'].split(':')[1]
                        exons.append(Exon(
                            transcript=transcripts[parent],
                            start=record['start'],
                            end=record['end'],
                        ))
                    if record['type'] == 'CDS':
                        parent = attributes['Parent'].split(':')[1]

                        if parent in genes and parent not in transcripts:
                            transcripts[parent] = Transcript.objects.create(
                                gene=genes[attributes['Parent']].split(':')[1],
                                transcript_id=attributes['Parent'].split(':')[1],
                                name=name,
                                start=record['start'],
                                end=record['end'],
                                transcript_type=genes[parent].gene_type,
                                source=record['source']
                            )
                        coding_sequences.append(CodingSequence(
                            transcript=transcripts[parent],
                            start=record['start'],
                            end=record['end'],
                        ))
                myprint('insert exons')
                Exon.objects.bulk_create(exons, batch_size=batch_size)
                myprint('insert coding_sequences')
                CodingSequence.objects.bulk_create(coding_sequences, batch_size=batch_size)

    def extract_sequence(self, chromosome:str, start: int, end: int) -> typing.Union[str, str]:
        return extract_sequence(self.assembly, chromosome, '+', start, end)


class SequenceObjectMixin:
    @property
    def sequence(self):
        return self.extract_sequence()

    @property
    def upstream(self):
        return self.extract_upstream()

    @property
    def downstream(self):
        return self.extract_downstream()

    def extract_upstream(self, length=SEQUENCE_PADDING):
        return self.extract_sequence(length)[:length]

    def extract_downstream(self, length=SEQUENCE_PADDING):
        return self.extract_sequence(length)[-length:]


class Gene(models.Model, SequenceObjectMixin):
    """Genes within design are represented by this model.

    Genes are stretches of DNA sequences that can be selected for editing.
    Genes can also contain one or multiple Transcripts, that can be selected for editing.

    organism -- foreign key to Organism table
    gene_id -- unique identifier for the gene, eg. Ensemb id
    name -- name of the gene eg. a HGNC symbol
    chromosome -- the scaffold where the sequence is located
    strand -- + or -
    start -- position of the first nucleotide
    end -- position of the last nucleotide
    gene_type -- type of the gene eg. protein_coding
    source -- which kind of curation was used to generate this gene annotation.
    """
    organism = models.ForeignKey('Organism', on_delete=models.CASCADE, related_name='genes')
    gene_id = models.CharField(max_length=200, db_index=True)
    name = models.CharField(max_length=200, db_index=True)
    chromosome = models.CharField(max_length=200)
    strand = models.CharField(max_length=1)
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()
    gene_type = models.CharField(max_length=200)
    source = models.CharField(max_length=200)

    class Meta:
        app_label = 'design'
        ordering = ('id', )

    def __str__(self):
        return f'{self.name} - {self.gene_id}'

    @property
    def degenerate_sequence(self):
        return []

    def extract_sequence(self, padding=0):
        return extract_sequence(self.organism.assembly, self.chromosome, self.strand, self.start-padding-1, self.end+padding)[0]

    def contextualize_alteration(self, position):
        """Wrapper for edit"""
        return position + SEQUENCE_PADDING


class Transcript(models.Model, SequenceObjectMixin):
    """Transcripts within design are represented by this model.

    Transcripts are stretches of DNA that can be edited.
    Transcripts can have a coding sequences, in which case edits can be specified by desired amino acid changes.

    gene -- the gene for which the transcript belongs
    transcript_id -- unique identifier for the transcript, eg. ensembl id
    name -- name of the transcript
    start -- position of the first nucleotide
    end -- position of the last nucleotide
    transcript_type -- type of transcript eg. protein_coding
    source -- how the transcript annotation was curated
    """
    gene = models.ForeignKey('Gene', on_delete=models.CASCADE, related_name='transcripts')
    transcript_id = models.CharField(max_length=200, db_index=True)
    name = models.CharField(max_length=200, db_index=True)
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()
    transcript_type = models.CharField(max_length=200)
    source = models.CharField(max_length=200)

    class Meta:
        app_label = 'design'
        ordering = ('transcript_id', )

    @property
    def organism(self):
        return self.gene.organism

    @property
    def codon_usage_table(self):
        return load_codon_table(self.gene.organism.codon_table)

    @property
    def strand(self):
        return self.gene.strand

    @property
    def coding_sequence(self):
        order_by = 'start'
        if self.strand == '-':
            order_by = '-start'
        return ''.join(cds.extract_sequence() for cds in self.coding_sequences.order_by(order_by).all())

    @property
    def translation(self):
        return translate(self.coding_sequence)

    def extract_sequence(self, padding=0):
        return extract_sequence(self.gene.organism.assembly, self.gene.chromosome, self.gene.strand, self.start-padding-1, self.end+padding)[0]

    def cds_to_transcript_position(self, cds_position):
        distance = 0
        order_by = 'start'
        if self.strand == '-':
            order_by = '-start'
        for cds in self.coding_sequences.order_by(order_by).all():
            length = cds.end - cds.start + 1
            if cds_position < distance + length:
                if self.gene.strand == '+':
                    start = cds.start - self.start
                else:
                    start = self.end - cds.end
                return start + cds_position - distance
            distance += length
        return None

    @staticmethod
    def cds_to_protein_position(cds_position):
        return cds_position // 3

    @staticmethod
    def protein_to_cds_position(protein_position):
        return protein_position*3, protein_position*3+2

    def protein_to_transcript_position(self, protein_position, insertions_before=0):
        start, end = self.protein_to_cds_position(protein_position-insertions_before)
        positions = [self.cds_to_transcript_position(position) for position in list(range(start, end+1))]
        return [p+insertions_before*3 for p in positions]

    @property
    def degenerate_sequence(self):
        dgn = ''.join([DEGENERATE_CODON_TABLE[aa] for aa in self.translation])
        offset = 0
        l = []
        order_by = 'start'
        if self.strand == '-':
            order_by = '-start'
        for cds in self.coding_sequences.order_by(order_by).all():
            length = cds.end - cds.start
            l.append((dgn[offset:offset + length + 1], self.cds_to_transcript_position(offset)))
            offset += length + 1
        return l

    def contextualize_alteration(self, position):
        return position + SEQUENCE_PADDING


class Exon(models.Model):
    transcript = models.ForeignKey('Transcript', on_delete=models.CASCADE, related_name='exons')
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()

    class Meta:
        app_label = 'design'

    def __str__(self):
        return f'{str(self.transcript)} Exon {self.start}-{self.end}'

    def extract_sequence(self):
        return extract_sequence(self.transcript.gene.organism.assembly, self.transcript.gene.chromosome, self.transcript.gene.strand, self.start-1, self.end)[0]


class CodingSequence(models.Model):
    transcript = models.ForeignKey('Transcript', on_delete=models.CASCADE, related_name='coding_sequences')
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()

    class Meta:
        app_label = 'design'

    def extract_sequence(self):
        return extract_sequence(self.transcript.gene.organism.assembly, self.transcript.gene.chromosome, self.transcript.gene.strand, self.start-1, self.end)[0]
