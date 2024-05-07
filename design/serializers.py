from Bio import Entrez
import django
from rest_framework import serializers

from .edits import EDITS
from .oligos import SequenceObject
from .models import Organism, Gene, Exon, Transcript, CodingSequence, SEQUENCE_PADDING

Entrez.email = django.conf.settings.DESIGN_ENTREZ_EMAIL


class OrganismSerializer(serializers.ModelSerializer):

    scaffolds = serializers.SerializerMethodField()

    class Meta:
        model = Organism
        fields = ('id', 'name', 'annotation', 'assembly', 'source', 'sequence_search', 'codon_table', 'scaffolds')

    def get_scaffolds(self, obj):
        return obj.scaffolds.split()


class ExonSerializer(serializers.ModelSerializer):

    class Meta:
        model = Exon
        fields = ('start', 'end')


class CodingSequenceSerializer(serializers.ModelSerializer):

    class Meta:
        model = CodingSequence
        fields = ('start', 'end')


class GeneTranscriptsSerializer(serializers.HyperlinkedModelSerializer):

    coding_sequences = serializers.SerializerMethodField()
    exons = serializers.SerializerMethodField()

    class Meta:
        model = Transcript
        fields = ('id', 'name', 'start', 'end', 'transcript_id', 'transcript_type', 'exons', 'coding_sequences')

    def get_exons(self, obj):
        return ExonSerializer(sorted([ex for ex in obj.exons.all()], key=lambda ex: ex.start,
                                     reverse=obj.strand == '-'), many=True).data

    def get_coding_sequences(self, obj):
        return CodingSequenceSerializer(sorted([cds for cds in obj.coding_sequences.all()], key=lambda cds: cds.start,
                                               reverse=obj.strand == '-'), many=True).data


class TranscriptSerializer(GeneTranscriptsSerializer):

    sequence = serializers.CharField(source='extract_sequence')
    organism_id = serializers.IntegerField(source='gene.organism_id')

    class Meta:
        model = Transcript
        fields = ('id', 'name', 'start', 'end', 'transcript_id', 'transcript_type', 'exons', 'coding_sequences',
                  'sequence', 'upstream', 'downstream', 'organism_id', 'strand', 'source', 'degenerate_sequence', 'translation')


class GeneSerializer(serializers.HyperlinkedModelSerializer):

    transcripts = GeneTranscriptsSerializer(many=True)
    sequence = serializers.CharField(source='extract_sequence')

    class Meta:
        model = Gene
        fields = ('id', 'organism_id', 'gene_id', 'name', 'chromosome', 'strand', 'start', 'end', 'sequence', 'upstream', 'downstream', 'transcripts',
                  'gene_type', 'source', 'url')


class GeneListSerializer(serializers.HyperlinkedModelSerializer):

    class Meta:
        model = Gene
        fields = ('id', 'gene_id', 'name', 'url')


class EditSerializer(serializers.Serializer):
    name = serializers.CharField(source='__name__')
    allow_genomic = serializers.BooleanField()
    default_options = serializers.SerializerMethodField()
    docstring = serializers.SerializerMethodField()
    required = serializers.ListField(
        child=serializers.CharField()
    )
    options = serializers.SerializerMethodField()

    @staticmethod
    def get_default_options(obj):
        return {k: (v[0].__name__, v[1]) for k, v in obj.default_options.items()}

    def get_options(self, obj):
        return {k: {
            'type': self.type2html(v),
            'value': v[1],
            'pattern': v[2],
            'required': k in obj.required
        } for k, v in obj.default_options.items()}

    @staticmethod
    def get_docstring(obj):
        docstring = obj.__doc__.splitlines()
        indent = min([len(line)-len(line.lstrip()) for line in docstring[1:] if line.lstrip()])
        docstring = '\n'.join([docstring[0]] + [line[indent:] for line in docstring[1:]])
        return docstring

    @staticmethod
    def type2html(typ):
        if type(typ[2]) is list:
            return 'select'
        if typ[0].__name__ == 'str':
            return 'text'
        if typ[0].__name__ == 'int':
            return 'number'
        return typ

    @staticmethod
    def parse_edit_dict(edit_dict: dict, organism: int, options: dict = None):
        if options is None:
            options = {}
        try:
            edit = EDITS[edit_dict['edit']]
        except KeyError:
            raise ValueError(f'Unknown edit type {edit_dict["edit"]}')
        sequence_type = edit_dict['sequence_type']
        sequence = edit_dict['sequence']
        organism = Organism.objects.get(pk=organism)
        option_string = edit_dict['options']

        padding = 0
        if sequence_type == 'transcript':
            sequence_objects = Transcript.objects.filter(name=sequence, gene__organism=organism).all()
            padding = SEQUENCE_PADDING
        elif sequence_type == 'gene':
            sequence_objects = Gene.objects.filter(name=sequence, organism=organism).all()
            padding = SEQUENCE_PADDING
        elif sequence_type == 'genomic':
            if not edit.allow_genomic:
                raise ValueError(f'edit {edit.__name__} is not allowed on genomic sequences')
            parsed_options = edit.parse_option_string(option_string)
            position = int(parsed_options['position']) - 1
            start = 0
            if position > 1000:
                parsed_options['position'] = 1001
                padding = 1001 - position
                start = position - 1000
            option_string = ','.join([f'{k}={e}' for k, e in parsed_options.items()])
            ref = organism.extract_sequence(sequence, start, position + 1000)[0]
            sequence_objects = [SequenceObject(ref, organism, name=sequence)]
        elif sequence_type == 'custom':
            if not edit.allow_genomic:
                raise ValueError(f'edit {edit.__name__} is not allowed on custom sequences')
            try:
                name, sequence = sequence.split(',')
            except ValueError:
                name = 'custom_sequence'
            sequence_objects = [SequenceObject(sequence, organism, name=name)]
        else:
            raise Exception(
                f"Unknown sequence_type {sequence_type}, accepted values: 'transcript, gene, genomic, custom'")
        if len(sequence_objects) > 1:
            raise Exception(f"Unambigious sequence identifier {sequence}, for organism {organism}")
        if len(sequence_objects) < 1:
            raise Exception(f"Unknown sequence identifier {sequence}, for organism {organism}")
        sequence_object = sequence_objects[0]
        edit = edit(sequence_object, option_string, **options)
        return edit, sequence_object, padding

    @classmethod
    def _validate(cls, edit_dict: dict) -> None:
        options = edit_dict['options']
        organism = edit_dict.pop('organism')
        edit, sequence_objects, padding = cls.parse_edit_dict(edit_dict, organism)
        if edit_dict['sequence_type'] == 'genomic':
            parsed_options = edit.parse_option_string(options)
            parsed_options['position'] = 1001
            options = ','.join([f'{k}={e}' for k, e in parsed_options.items()])
        options = edit.parse_option_string(options)
        edit.validate_options(options, sequence_objects)

    @classmethod
    def validate_input(cls, edit_dict: dict) -> dict:
        try:
            cls._validate(edit_dict)
        except ValueError as e:
            raise serializers.ValidationError(str(e))
        return edit_dict

    @classmethod
    def validate_input_list(cls, edits_dict: dict):
        organism = edits_dict.pop('organism')
        edits = []
        errors = []
        for i, edit in enumerate(edits_dict['edits'], 1):
            edit.update({'organism': organism})
            try:
                cls._validate(edit)
            except ValueError as e:
                edit.update({'error': str(e), 'line': i})
                errors.append(edit)
            else:
                if edit in edits:
                    edit.update({'error': 'Edit already added', 'line': i})
                    errors.append(edit)
                else:
                    edits.append(edit)
        return edits, errors

    @classmethod
    def validate_clinvar_list(cls, clinvar_list: dict):
        ids = [e['clinvar_id'].replace('VCV', '') for e in clinvar_list]
        repair = [e['repair'] for e in clinvar_list]
        edits = []
        errors = []
        handle = Entrez.esummary(db="clinvar", id=','.join(ids), retmode='xml')
        records = Entrez.read(handle, validate=False)
        for i, record in enumerate(records['DocumentSummarySet']['DocumentSummary']):
            edit = {'clinvar_id': f'VCV{ids[i]}', 'repair': repair[i]}
            variation_set = record['variation_set'][0]
            try:
                chromosome, position, from_, to_ = variation_set['canonical_spdi'].split(':')
            except (ValueError, KeyError):
                edit.update({'error': 'Could not find/parse Canonical SPDI', 'line': i})
                errors.append(edit)
            else:
                if edit in edits:
                    edit.update({'error': 'Edit already added', 'line': i})
                    errors.append(edit)
                else:
                    edits.append(edit)

        return edits, errors


class NucleaseSerializer(serializers.Serializer):
    name = serializers.CharField(source='__name__')
    pam = serializers.CharField(source='pam_motif')
    spacer_length = serializers.IntegerField()
    cut_site_position = serializers.IntegerField(source='_cut_site_position')
    scaffolds = serializers.DictField(child=serializers.CharField())
    cloning_strategies = serializers.SerializerMethodField()
    can_score_spacers = serializers.BooleanField()
    docstring = serializers.SerializerMethodField()
    options = serializers.SerializerMethodField()

    @staticmethod
    def remove_indent(string):
        string = string.splitlines()
        if len(string) > 1:
            indent = min([len(line) - len(line.lstrip()) for line in string[1:] if line.lstrip()])
            return '\n'.join([string[0]] + [line[indent:] for line in string[1:]])
        return string[0]

    def get_docstring(self, obj):
        return self.remove_indent(obj.__doc__)

    def get_cloning_strategies(self, obj):
        return [(k,
                 self.remove_indent(v.help_text(obj.scaffolds[obj.default_scaffold])),
                 {k2: self.opt2dict(v2) for k2, v2 in v.options.items()}) for k, v in obj.cloning_strategies.items()]

    @staticmethod
    def type2html(typ):
        if type(typ[2]) is list:
            return 'select'
        if typ[0].__name__ == 'str':
            return 'text'
        if typ[0].__name__ == 'int':
            return 'number'
        if typ[0].__name__ == 'bool':
            return 'checkbox'
        return typ

    def opt2dict(self, opt):
        return {
            'type': self.type2html(opt),
            'value': opt[1],
            'pattern': opt[2],
            'help': opt[3],
        }

    def get_options(self, obj):
        return {k: self.opt2dict(v) for k, v in obj.options.items()}
