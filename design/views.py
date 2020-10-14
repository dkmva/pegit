import base64
import json
import os

from Bio import Entrez
import django.conf
from django.http import FileResponse
from django.shortcuts import render
from rest_framework import viewsets, filters, status, serializers
from rest_framework.decorators import action
from rest_framework.exceptions import NotFound
from rest_framework.response import Response

from prime.celery import app as celery_app
from design import EDITS, NUCLEASES
from design.serializers import EditSerializer, GeneSerializer, GeneListSerializer, NucleaseSerializer,\
    OrganismSerializer, TranscriptSerializer
from design.models import Organism, Gene, Transcript
from design.interface import Job, JobSerializer, create_oligos_chain

Entrez.email = django.conf.settings.DESIGN_ENTREZ_EMAIL


def index(request):
    return render(request, 'design/index.html')


class MultiSerializerViewSet(viewsets.ReadOnlyModelViewSet):
    serializers = {
        'default': None,
    }

    def get_serializer_class(self):
        return self.serializers.get(self.action, self.serializers['default'])


class OrganismViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Organism.objects.all()
    serializer_class = OrganismSerializer

    @action(methods=['GET'], url_path='region/(?P<region>.+:\d+-\d+)', detail=True)
    def region(self, request, pk, region, *args, **kwargs):
        organism = self.get_object()

        chromosome, start, end = region.replace(':', '-').split('-')

        sequence, error = organism.extract_sequence(chromosome, int(start)-1, int(end))
        if error:
            return Response({'message': error}, status.HTTP_404_NOT_FOUND)

        return Response({'sequence': sequence})


class GeneViewSet(MultiSerializerViewSet):
    pagination_class = None
    serializers = {
        'default': GeneSerializer,
        'list': GeneListSerializer,
    }

    filter_backends = (filters.SearchFilter, )
    search_fields = ['gene_id', '^name']

    def get_queryset(self):
        queryset = Gene.objects.select_related('organism').prefetch_related('transcripts',
                                                                            'transcripts__coding_sequences',
                                                                            'transcripts__exons').all()
        if 'organism' in self.request.query_params:
            queryset = queryset.filter(organism__id=self.request.query_params['organism'])

        return queryset


class TranscriptViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Transcript.objects.select_related('gene').prefetch_related('exons', 'coding_sequences').all()
    serializer_class = TranscriptSerializer


class EditViewSet(viewsets.ViewSet):

    serializer_class = EditSerializer

    def list(self, request, *args, **kwargs):
        serializer = EditSerializer(instance=EDITS.values(), many=True)
        return Response(serializer.data)

    def retrieve(self, request, pk, *args, **kwargs):
        serializer = EditSerializer(instance=EDITS[pk])
        return Response(serializer.data)

    @action(methods=['POST'], detail=False)
    def validate(self, request, *args, **kwargs):
        try:
            EditSerializer.validate_input(request.data)
        except serializers.ValidationError as e:
            return Response(serializers.as_serializer_error(e), status=status.HTTP_400_BAD_REQUEST)
        return Response({})

    @action(methods=['POST'], detail=False)
    def validate_list(self, request, *args, **kwargs):
        try:
            edits, errors = EditSerializer.validate_input_list(request.data)
        except Exception as e:
            return Response({'error': str(e)}, status=status.HTTP_400_BAD_REQUEST)
        return Response({'edits': edits, 'errors': errors})

    @action(methods=['GET'], detail=False)
    def advanced_options(self, request, *args, **kwargs):
        return Response(django.conf.settings.DESIGN_CONF['default_options'])


class NucleaseViewSet(viewsets.ViewSet):

    serializer_class = NucleaseSerializer

    def list(self, request, *args, **kwargs):
        serializer = NucleaseSerializer(instance=NUCLEASES.values(), many=True)
        return Response(serializer.data)


class JobViewSet(viewsets.ViewSet):

    def list(self, request, *args, **kwargs):
        # Jobs should not be viewable by everyone
        # In the future, consider restricting access and provide thin serializer
        return Response([])

    def retrieve(self, request, pk, *args, **kwargs):
        try:
            return FileResponse(open(os.path.join(django.conf.settings.DESIGN_OUTPUT_FOLDER, pk, 'summary.json'), 'rb'))
        except FileNotFoundError:
            raise NotFound

    @action(detail=True, methods=['GET'], url_path='edit(?P<edit>[a-z0-9]+)')
    def edit(self, request, pk, edit, *args, **kwargs):
        try:
            return FileResponse(open(os.path.join(django.conf.settings.DESIGN_OUTPUT_FOLDER, pk, f'edit{edit}.json'), 'rb'))
        except FileNotFoundError:
            raise NotFound

    @action(detail=True, methods=['GET'])
    def download(self, request, pk, *args, **kwargs):
        jobdir = os.path.join(django.conf.settings.DESIGN_OUTPUT_FOLDER, pk)
        xlsxfile = [e for e in os.listdir(jobdir) if e.endswith('xlsx')][0]
        return FileResponse(open(os.path.join(jobdir, xlsxfile), 'rb'))

    @action(detail=True, methods=['GET'])
    def queue_position(self, request, pk, *args, **kwargs):
        with celery_app.pool.acquire(block=True) as conn:
            tasks = conn.default_channel.client.lrange('design_queue', 0, -1)
        pks = [json.loads(base64.b64decode(json.loads(t)['body']))[0][0] for t in tasks][::-1]

        try:
            return Response({'position': pks.index(pk) + 1})
        except ValueError:
            return Response({'position': ''})

    def create(self, request, *args, **kwargs):
        pk = request.data.get('organism', None)
        edits = request.data.get('edits', None)
        advanced_options = request.data.get('advanced_options', None)
        nuclease = request.data.get('nuclease', None)
        organism = Organism.objects.get(pk=pk)
        j = Job(organism, options=advanced_options, edits=edits, nuclease=nuclease)
        j.save()

        create_oligos_chain(j.job_id)

        return Response(JobSerializer(j).data, status=status.HTTP_201_CREATED)

    @action(detail=False, methods=['POST'])
    def clinvar(self, request, *args, **kwargs):
        edits = request.data.get('edits', None)

        organism = Organism.objects.get(assembly=django.conf.settings.DESIGN_CLINVAR_ORGANISM)
        advanced_options = request.data.get('advanced_options', None)
        nuclease = request.data.get('nuclease', None)
        j = Job(organism, options=advanced_options, nuclease=nuclease)
        j.edits = j.clinvar2edit(edits)
        j.save()

        create_oligos_chain(j.job_id)

        return Response(JobSerializer(j).data, status=status.HTTP_201_CREATED)


class ClinvarSearch(viewsets.ViewSet):

    def create(self, request, *args, **kwargs):
        query = request.data.get('query', None)
        if query:
            with Entrez.esearch(db="clinvar", retmax=30, term=query, retmode='xml') as handle:
                ids = Entrez.read(handle)['IdList']
            with Entrez.esummary(db="clinvar", id=','.join(ids), retmode='xml') as handle:
                records = Entrez.read(handle)['DocumentSummarySet']['DocumentSummary']

            return Response({'results': records})
        return Response({'results': []})

    @action(methods=['POST'], detail=False)
    def validate_list(self, request, *args, **kwargs):
        edits, errors = EditSerializer.validate_clinvar_list(request.data)
        return Response({'edits': edits, 'errors': errors})
