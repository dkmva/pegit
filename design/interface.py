from __future__ import annotations
import copy
from collections import defaultdict
import dataclasses
import itertools
import enum
import json
import os
import typing
import uuid

from Bio import Entrez
from celery import shared_task
from celery.exceptions import SoftTimeLimitExceeded
import django
from djangorestframework_camel_case.util import camelize, underscoreize
from rest_framework import serializers
import pandas as pd

from design.helpers import check_primer_specificity
from design.models import Organism
from design.nucleases import NUCLEASES
from design.serializers import OrganismSerializer, EditSerializer

Entrez.email = django.conf.settings.DESIGN_ENTREZ_EMAIL


class JobStatus(enum.Enum):
    QUEUED = 'Queued'
    FINDING_PEGRNAS = 'Finding pegRNAs'
    QUEUED_SPECIFITY = 'Queued for specificity checks'
    CHECKING_PRIMER_SPECIFICITY = 'Checking primer specificity'
    CHECKING_SGRNA_SPECIFICITY = 'Checking sgRNA specificity'
    COMPLETED_STATUS = 'Completed'
    FAILED_STATUS = 'Failed'

    def __repr__(self):
        return self.value

    def __str__(self):
        return self.value


def make_uuid_string() -> str:
    return str(uuid.uuid4())


@dataclasses.dataclass
class Job:
    organism: Organism
    options: dict = dataclasses.field(default_factory=dict)
    job_name: str = 'pegIT'
    nuclease: str = dataclasses.field(default=None)
    status: JobStatus = JobStatus.QUEUED
    job_id: str = dataclasses.field(default_factory=make_uuid_string)
    edits: list = dataclasses.field(default_factory=list)
    summary: list = dataclasses.field(default_factory=list)
    warning: str = dataclasses.field(default=None)
    jobdir: str = dataclasses.field(default=None, init=False)
    output_folder: dataclasses.InitVar(typing.Optional[str]) = None

    def __post_init__(self, output_folder: typing.Optional[str]) -> None:
        options = django.conf.settings.DESIGN_CONF['default_options'].copy()
        options.update(self.options)
        self.options = options

        if self.nuclease is None:
            self.nuclease = django.conf.settings.DESIGN_CONF['default_nuclease']

        if output_folder is None:
            output_folder = django.conf.settings.DESIGN_OUTPUT_FOLDER
        self.jobdir = os.path.join(output_folder, self.job_id)

    @property
    def results(self):
        for i, edit in enumerate(self.edits):
            try:
                with open(os.path.join(self.jobdir, f'edit{i}.json')) as f:
                    data = json.load(f)
                    primers = data.pop('primers')
                    data = underscoreize(data)
                    data['pegRNAs'] = data.pop('peg_rn_as')
                    data['primers'] = primers
                    yield data
            except FileNotFoundError:
                break

    def save_result(self, result, edit_index):
        with open(os.path.join(self.jobdir, f'edit{edit_index}.json'), 'w') as f:
            json.dump(camelize(result), f, indent=4)

    def save(self, as_json=True, as_excel=True, edit_index=None) -> None:
        try:
            os.makedirs(self.jobdir)
        except FileExistsError:
            pass
        if as_json:
            if edit_index is None or edit_index % 10 == 0:
                with open(os.path.join(self.jobdir, 'summary.json'), 'w') as f:
                    d = camelize(JobSerializer(self).data)
                    del d['organism']['scaffolds']
                    json.dump(d, f, indent=4)
        if as_excel:
            self.export_excel()

    @classmethod
    def load_from_disk(cls, job_id: str) -> Job:
        jobdir = os.path.join(django.conf.settings.DESIGN_OUTPUT_FOLDER, job_id)
        with open(os.path.join(jobdir, 'summary.json')) as f:
            data = underscoreize(json.load(f), ignore_fields=['pegRNAs'])
            data['organism'] = Organism.objects.get(pk=data['organism']['id'])
            for s in data['summary']:
                s['pegRNAs'] = s.pop('peg_rn_as')
            job = cls(**data, output_folder=django.conf.settings.DESIGN_OUTPUT_FOLDER)
            for i, edit in enumerate(job.edits):
                try:
                    with open(os.path.join(jobdir, f'edit{i}.json')) as f:
                        data = underscoreize(json.load(f))
                        data['pegRNAs'] = data.pop('peg_rn_as')
                except FileNotFoundError:
                    break
            return job

    def import_edit_list(self, edit_list: str) -> None:
        """Import edits from a TSV file"""
        header = ['sequence_type', 'sequence', 'edit', 'options']
        self.edits = []
        with open(edit_list, 'r') as f:
            for line in f:
                self.edits.append(dict(zip(header, line.split())))

    def save_edit_list(self, filepath: str) -> None:
        """Save edits to a TSV file"""
        with open(filepath, 'w') as f:
            for edit in self.edits:
                sequence_type = edit['sequence_type']
                sequence = edit['sequence']
                options = edit['options']
                edit = edit['edit']
                f.write(f'{sequence_type}\t{sequence}\t{edit}\t{options}\n')

    def import_clinvar_list(self, edit_list: str) -> None:
        """Import edits from a clinVar TSV file"""
        if self.organism.assembly != django.conf.settings.DESIGN_CLINVAR_ORGANISM:
            raise ValueError('Organism does not correspond to a ClinVar DB')
        header = ['clinvar_id', 'repair']
        clinvar_edits = []
        with open(edit_list, 'r') as f:
            for line in f:
                clinvar_edits.append(dict(zip(header, line.split())))
        self.edits = self.clinvar2edit(clinvar_edits)

    def save_clinvar_list(self, filepath: str) -> None:
        """Save edits to a clinvar TSV file"""
        if self.organism.assembly != django.conf.settings.DESIGN_CLINVAR_ORGANISM:
            raise ValueError('Organism does not correspond to a ClinVar DB')
        with open(filepath, 'w') as f:
            for edit in self.edits:
                clinvar_id = edit['sequence']
                repair = edit['options'].split('=')[1]
                f.write(f'{clinvar_id}\t{repair}\n')

    @staticmethod
    def clinvar2edit(clinvar_list: typing.List[dict]) -> list:
        """Convert a list of clinvar edits to regular edits"""
        ids = [e['clinvar_id'].replace('VCV', '') for e in clinvar_list]
        repair = [e['repair'] for e in clinvar_list]

        handle = Entrez.esummary(db="clinvar", id=','.join(ids), retmode='xml')
        records = Entrez.read(handle)
        edits = []
        for i, record in enumerate(records['DocumentSummarySet']['DocumentSummary']):
            variation_set = record['variation_set'][0]
            chromosome, position, from_, to_ = variation_set['canonical_spdi'].split(':')

            edit_dict = {
                'sequence_type': 'genomic',
                'sequence': chromosome,
                'edit': 'Substitution',
                'options': f'from={from_},to={to_},position={int(position)+1}'
            }
            if repair[i] == 'true':
                edit_dict['options'] += ',repair=true'

            edits.append(edit_dict)
        return edits

    def export_excel(self) -> None:
        """Export jobdata to excel"""
        writer = pd.ExcelWriter(os.path.join(self.jobdir, f'{self.job_name}.xlsx'))
        wb = writer.book

        # SUMMARY
        summary = pd.DataFrame(self.summary)
        summary.to_excel(writer, 'Summary', index=False, startrow=8)
        ws = wb['Summary']
        ws[f'A8'] = 'Edits'
        ws[f'A8'].style = 'Headline 1'

        organism = OrganismSerializer(self.organism).data
        del organism['scaffolds']

        pd.DataFrame([organism]).to_excel(writer, 'Summary', index=False, startrow=1)
        ws[f'A1'] = 'Organism'
        ws[f'A1'].style = 'Headline 1'

        pd.DataFrame([self.options]).to_excel(writer, 'Summary', index=False, startrow=5)
        ws[f'A5'] = 'Options'
        ws[f'A5'].style = 'Headline 1'

        writer.save()

        summary_pegrnas = []
        summary_primers = []

        pegRNA_header = ['#pegRNA', 'spacer', 'score', 'pam_disrupted', 'pam_silenced', 'distance', 'strand',
                         'nuclease', 'pbs_length', 'rt_template_length', 'pbs', 'rt_template', 'offtargets',
                         'spacer_oligo_top', 'spacer_oligo_bottom', 'extension_oligo_top', 'extension_oligo_bottom']
        # EDITS
        for i, result in enumerate(self.results, 1):
            wsname = f"{i} {self.edits[i-1]['sequence']} {self.edits[i-1]['edit']}"
            pegRNAs = copy.deepcopy(result['pegRNAs'])
            if not pegRNAs:
                summary_pegrnas.append({"#Edit": i})
                summary_primers.append({"#Edit": i})
                continue
            all_nicking = []
            all_alternate = []
            for j, p in enumerate(pegRNAs):
                p['#pegRNA'] = j+1
                try:
                    p['offtargets'] = ','.join(str(e) for e in p['offtargets'][0])
                except KeyError:
                    p['offtargets'] = 'unknown'
                oligos = p.pop('oligos')
                p['spacer_oligo_top'] = oligos['spacer']['top']
                p['spacer_oligo_bottom'] = oligos['spacer']['bottom']
                p['extension_oligo_top'] = oligos['extension']['top']
                p['extension_oligo_bottom'] = oligos['extension']['bottom']
                nicking = p.pop('nicking')
                alternate_extensions = p.pop('alternate_extensions')
                pegRNAs[j] = {k: p[k] for k in pegRNA_header}
                try:
                    pegRNAs[j]['nsgRNA_oligo_top'] = nicking[0]['top']
                    pegRNAs[j]['nsgRNA_oligo_bottom'] = nicking[0]['bottom']
                except IndexError:
                    pegRNAs[j]['nsgRNA_oligo_top'] = ''
                    pegRNAs[j]['nsgRNA_oligo_bottom'] = ''
                for n in nicking:
                    n['#pegRNA'] = j + 1
                    try:
                        n['offtargets'] = ','.join(str(e) for e in n['offtargets'][0])
                    except KeyError:
                        n['offtargets'] = 'unknown'
                    n['oligo_top'] = n['top']
                    n['oligo_bottom'] = n['bottom']
                    n = {k: n[k] for k in ['#pegRNA', 'kind', 'position', 'spacer', 'score', 'wt_score', 'offtargets',
                                             'oligo_top', 'oligo_bottom']}
                    all_nicking.append(n)

                for ext in alternate_extensions:
                    oligos = ext.pop('oligos')
                    ext['#pegRNA'] = j + 1
                    ext['oligo_top'] = oligos['top']
                    ext['oligo_bottom'] = oligos['bottom']
                    ext = {k: ext[k] for k in ['#pegRNA', 'pbs_length', 'rt_template_length', 'sequence', 'pbs_gc',
                                               'rt_gc', 'oligo_top', 'oligo_bottom',]}
                    all_alternate.append(ext)

            primers = [p['primers'] for p in result['primers']]
            all_primers = []
            for j, p in enumerate(primers):
                d = {}
                for k in p.keys():
                    for kk in p[k]:
                        d[f'{k}_{kk}'] = p[k][kk]
                all_primers.append(d)
                if j == 0:
                    summary_dict = {'#Edit': i, 'left_sequence': '', 'right_sequence': ''}
                    summary_dict.update({k: d[k] for k in sorted(d, key=lambda x: (x.startswith('p'), x.startswith('r')))})
                    summary_primers.append(summary_dict)

            start_row = 1

            for j, df in enumerate([pegRNAs, all_nicking, all_primers, all_alternate]):
                df = pd.DataFrame(df)
                df.to_excel(writer, wsname, index=False, startrow=start_row)
                ws = wb[wsname]
                ws[f'A{start_row}'] = ['pegRNAs', 'nsgRNAs', 'primers', 'alternate extensions'][j]
                ws[f'A{start_row}'].style = 'Headline 1'
                start_row += len(df.index) + 3
            summary_dict = {'#Edit': i}
            summary_dict.update(pegRNAs[0])
            del summary_dict["#pegRNA"]
            summary_pegrnas.append(summary_dict)

        start_row = 8
        ws = wb['Summary']
        start_row += len(summary.index) + 3
        pd.DataFrame(summary_pegrnas).to_excel(writer, 'Summary', index=False, startrow=start_row)
        ws[f'A{start_row}'].style = 'Headline 1'
        ws[f'A{start_row}'] = 'pegRNAs'
        start_row += len(summary_pegrnas) + 3
        pd.DataFrame(summary_primers).to_excel(writer, 'Summary', index=False, startrow=start_row)
        ws[f'A{start_row}'].style = 'Headline 1'
        ws[f'A{start_row}'] = 'Primers'
        writer.save()

    def create_oligos(self):
        """Design pegRNAs for job edits. Saves to jobdir folder"""
        try:
            self.status = JobStatus.FINDING_PEGRNAS
            self.summary = []
            self.save(as_excel=False)
            for i, edit in enumerate(self.edits):
                sequence_type = edit['sequence_type']
                sequence = edit['sequence']
                options = edit['options']

                edit_dict = edit.copy()
                edit_dict.update({'organism': self.organism.id})

                result = {'warning': None, 'sequence': sequence, 'sequence_type': sequence_type, 'options': options,
                          'pegRNAs': [], 'primers': [],
                          'sequence_object': {
                              'name': None,
                              'source': None,
                              'sequence': None
                          }}

                try:
                    edit, sequence_object, padding = EditSerializer.parse_edit_dict(edit_dict, self.organism.pk, dict({'nuclease': self.nuclease}, **self.options))
                    result.update(edit.create_oligos())

                    if len(result['pegRNAs']) < 1:
                        result['warning'] = 'No pegRNAs found'

                    result['start'] -= padding
                    result['edit'] = edit.__class__.__name__
                    result['sequence_object'] = SequenceObjectSerializer(sequence_object).data
                    if sequence_type == 'custom':
                        result['sequence'] = ','.join([sequence_object.name, sequence_object.sequence])
                except Exception as e:
                    result['warning'] = f'An unknown error occured, {repr(e)}'

                finally:
                    self.save_result(result, i)
                    self.summary.append({
                        'sequence': result['sequence'],
                        'sequence_type': sequence_type,
                        'edit': edit.__class__.__name__,
                        'options': options,
                        'pegRNAs': len(result['pegRNAs']),
                        'warning': result['warning']
                    })
                    self.save(as_excel=False, edit_index=i)

            self.status = JobStatus.CHECKING_PRIMER_SPECIFICITY
            self.save()

            all_primers = list(itertools.chain.from_iterable([result['primers']for result in self.results]))
            all_primers = check_primer_specificity(all_primers, self.organism.assembly, os.path.join(
                                                                           django.conf.settings.DESIGN_OUTPUT_FOLDER,
                                                                           self.job_id), **self.options)

            i = 0
            for j, result in enumerate(self.results):
                num_primers = len(result.get('primers', []))
                result['primers'] = sorted(all_primers[i:i+num_primers], key=lambda x: x['product_count'])
                i += num_primers
                self.save_result(result, j)

            self.status = JobStatus.CHECKING_SGRNA_SPECIFICITY
            self.save()

            spacers = defaultdict(lambda: dict())
            for result in self.results:
                for pegRNA in result['pegRNAs']:
                    spacers[pegRNA['nuclease']][pegRNA['spacer']] = 1
                    for nicking in pegRNA['nicking']:
                        spacers[pegRNA['nuclease']][nicking['spacer']] = 1

            for nuclease in spacers:
                spacers[nuclease] = list(spacers[nuclease].keys())
                counts, binders = NUCLEASES[nuclease].find_off_targets(spacers[nuclease],
                                                                       self.organism.assembly,
                                                                       os.path.join(
                                                                           django.conf.settings.DESIGN_OUTPUT_FOLDER,
                                                                           self.job_id))
                for i, result in enumerate(self.results):
                    for j, hits in enumerate(zip(counts, binders)):
                        for pegRNA in result['pegRNAs']:
                            if pegRNA['nuclease'] == nuclease and pegRNA['spacer'] == spacers[nuclease][j]:
                                pegRNA['offtargets'] = hits
                            for nicking in pegRNA['nicking']:
                                if pegRNA['nuclease'] == nuclease and nicking['spacer'] == spacers[nuclease][j]:
                                    nicking['offtargets'] = hits
                    self.save_result(result, i)

            self.status = JobStatus.COMPLETED_STATUS
            self.save()
        except Exception as e:
            self.status = JobStatus.FAILED_STATUS
            self.warning = f'An unknown error occured, {repr(e)}'
            self.save()

            return e


class JobSerializer(serializers.Serializer):
    organism = OrganismSerializer()
    status = serializers.CharField()
    job_name = serializers.CharField()
    job_id = serializers.CharField()
    summary = serializers.ListField(child=serializers.DictField(child=serializers.CharField()))
    edits = serializers.ListField(child=serializers.DictField(child=serializers.CharField()))
    options = serializers.DictField(child=serializers.IntegerField())
    nuclease = serializers.CharField(required=False)
    warning = serializers.CharField()


class SequenceObjectSerializer(serializers.Serializer):
    id = serializers.SerializerMethodField()
    name = serializers.CharField()
    source = serializers.CharField()

    def get_id(self, obj):
        try:
            return str(obj.transcript_id)
        except AttributeError:
            pass
        try:
            return str(obj.gene_id)
        except AttributeError:
            return None


@shared_task()
def create_oligos_background(job_id: str) -> None:
    """Used by celery to run create_oligos in background."""
    job = Job.load_from_disk(job_id)
    e = job.create_oligos()

    if isinstance(e, Exception):
        raise e
