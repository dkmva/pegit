from __future__ import annotations
import copy
from collections import defaultdict
import dataclasses
import itertools
import enum
import json
import os
import string
import typing
import uuid

from Bio import Entrez
import celery
import django
from djangorestframework_camel_case.util import camelize, underscoreize
from rest_framework import serializers
import pandas as pd

from design.helpers import check_primer_specificity
from design.models import Organism
from design.nucleases import NUCLEASES
from design.serializers import OrganismSerializer, EditSerializer

# Monkey patch pandas to write excel files by rows instead of columns.
# Enables writing excel files without putting everything into memory.
# https://github.com/pandas-dev/pandas/issues/34710
from pandas.io.formats.excel import ExcelFormatter, ExcelCell


DESIGN_EDIT_GROUP_SIZE = 5


def write_excel_by_rows(self, coloffset: int):
    if self.styler is None:
        styles = None
    else:
        styles = self.styler._compute().ctx
        if not styles:
            styles = None
    xlstyle = None
    for rowidx in range(self.df.shape[0]):
        for colidx in range(len(self.columns)):
            if styles is not None:
                xlstyle = self.style_converter(";".join(styles[rowidx, colidx]))
            yield ExcelCell(self.rowcounter + rowidx, colidx + coloffset, self.df.iloc[rowidx, colidx], xlstyle)


ExcelFormatter._generate_body = write_excel_by_rows
# Monkey patch done


Entrez.email = django.conf.settings.DESIGN_ENTREZ_EMAIL


class JobStatus(enum.Enum):
    QUEUED = 'Queued'
    FINDING_PEGRNAS = 'Finding pegRNAs'
    QUEUED_PRIMER_SPECIFITY = 'Queued for primer specificity'
    CHECKING_PRIMER_SPECIFICITY = 'Checking primer specificity'
    QUEUED_SGRNA_SPECIFITY = 'Queued for sgRNA specificity'
    CHECKING_SGRNA_SPECIFICITY = 'Checking sgRNA specificity'
    EXPORTING_EXCEL = 'Exporting Excel file'
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
    nuclease_options: dict = dataclasses.field(default_factory=dict)
    cloning_strategy: str = None
    cloning_options: dict = dataclasses.field(default_factory=dict)
    status: JobStatus = JobStatus.QUEUED
    job_id: str = dataclasses.field(default_factory=make_uuid_string)
    edits: list = dataclasses.field(default_factory=list)
    summary: list = dataclasses.field(default_factory=list)
    warning: typing.Optional[str] = dataclasses.field(default=None)
    jobdir: str = dataclasses.field(default=None, init=False)
    output_folder: dataclasses.InitVar(typing.Optional[str]) = None
    run_bowtie: bool = True
    design_primers: bool = True
    excel_exported: bool = False

    def __post_init__(self, output_folder: typing.Optional[str]) -> None:
        options = django.conf.settings.DESIGN_CONF['default_options'].copy()
        options.update(self.options)
        self.options = options

        if self.nuclease is None:
            self.nuclease = django.conf.settings.DESIGN_CONF['default_nuclease']

        if self.cloning_strategy is None:
            self.cloning_strategy = next(iter(NUCLEASES[self.nuclease].cloning_strategies))

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

    def save(self) -> None:
        os.makedirs(os.path.join(self.jobdir, 'excel_tmp'), exist_ok=True)
        with open(os.path.join(self.jobdir, 'summary.json'), 'w') as f:
            d = camelize(JobSerializer(self).data, ignore_fields=('cloning_options', 'nuclease_options'))
            del d['organism']['scaffolds']
            json.dump(d, f, indent=4)

    @classmethod
    def load_from_disk(cls, job_id: str) -> Job:
        jobdir = os.path.join(django.conf.settings.DESIGN_OUTPUT_FOLDER, job_id)
        with open(os.path.join(jobdir, 'summary.json')) as f:
            data = underscoreize(json.load(f), ignore_fields=('pegRNAs', 'cloning_options', 'nuclease_options'))
            data['organism'] = Organism.objects.get(pk=data['organism']['id'])
            for s in data['summary']:
                s['pegRNAs'] = s.pop('peg_rn_as')
            job = cls(**data, output_folder=django.conf.settings.DESIGN_OUTPUT_FOLDER)
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

    def clinvar2edit(self, clinvar_list: typing.List[dict], naming='accession') -> list:
        """Convert a list of clinvar edits to regular edits"""
        ids = [e['clinvar_id'].replace('VCV', '') for e in clinvar_list]
        repair = [e['repair'] for e in clinvar_list]

        handle = Entrez.esummary(db="clinvar", id=','.join(ids), retmode='xml')
        records = Entrez.read(handle)
        edits = []
        for i, record in enumerate(records['DocumentSummarySet']['DocumentSummary']):
            variation_set = record['variation_set'][0]
            try:
                chromosome, position, from_, to_ = variation_set['canonical_spdi'].split(':')
            except (ValueError, KeyError):
                print(f'Could not find/parse Canonical SPDI, line {i}')
                continue
            position = int(position) - 1
            parsed_position = position
            start = 0
            if position > 1000:
                parsed_position = 1001
                start = position - 1000

            ref = self.organism.extract_sequence(chromosome, start, position + 1000)[0]

            edit_dict = {
                'sequence_type': 'custom',
                'sequence': f'{record[naming]},{ref}',
                'edit': 'Substitution',
                'options': f'from={from_},to={to_},position={int(parsed_position)+1}'
            }
            if repair[i] == 'true':
                edit_dict['options'] += ',repair=true'

            edits.append(edit_dict)
        return edits

    def export_excel(self) -> None:
        """Export jobdata to excel"""

        cloning_strategy = NUCLEASES[self.nuclease].cloning_strategies[self.cloning_strategy]

        excel_name = ''.join([c if c in f"-_.() {string.ascii_letters}{string.digits}" else "_" for c in self.job_name])

        writer = pd.ExcelWriter(os.path.join(self.jobdir, f'{excel_name}.xlsx'), engine='xlsxwriter',
                                options=dict(constant_memory=True, tmpdir=os.path.join(self.jobdir, 'excel_tmp')))
        wb = writer.book
        heading = wb.add_format({'bold': True, 'font_size': 15})
        ws = wb.add_worksheet('Summary')
        writer.sheets['Summary'] = ws

        # SUMMARY
        summary_column_widths = defaultdict(lambda: 0)
        organism = OrganismSerializer(self.organism).data

        ws.write('A1', 'Organism', heading)
        ws.write('A2', 'Name')
        ws.write('B2', organism['name'])
        ws.write('A3', 'Source')
        ws.write('B3', organism['source'])

        summary_start_row = 5

        for title, opt in (('Options', self.options),
                    ('Nuclease options', self.nuclease_options),
                    ('Cloning options', self.cloning_options)):
            if opt:
                ws.write(f'A{summary_start_row}', title, heading)
                df = pd.DataFrame([self.options])
                df.to_excel(writer, 'Summary', index=False, startrow=summary_start_row)
                for col_idx, column in enumerate(df):
                    column_length = max(df[column].astype(str).map(len).max(), len(column))
                    summary_column_widths[col_idx] = max(column_length, summary_column_widths[col_idx])
                summary_start_row += 3

        summary = pd.DataFrame(self.summary)
        del summary['best']
        ws.write(f'A{summary_start_row}', 'Edits', heading)
        summary.to_excel(writer, 'Summary', index=False, startrow=summary_start_row)
        for col_idx, column in enumerate(summary):
            column_length = max(summary[column].astype(str).map(len).max(), len(column))
            summary_column_widths[col_idx] = max(column_length, summary_column_widths[col_idx])

        summary_pegrnas = []
        summary_primers = []

        pegRNA_header = ['#pegRNA', 'spacer', 'score', 'pam_disrupted', 'pam_silenced', 'distance', 'strand',
                         'nuclease', 'pbs_length', 'rt_template_length', 'pbs', 'rt_template', 'offtargets', 'pegRNA', 'nsgRNA']\
                        + cloning_strategy.excel_oligo_headers

        # EDITS
        for i, result in enumerate(self.results, 1):
            wsname = f"{i} {self.edits[i-1]['sequence'].split(',')[0]} {self.edits[i-1]['edit']}"[:31]
            pegRNAs = copy.deepcopy(result['pegRNAs'])
            if not pegRNAs:
                summary_pegrnas.append({"#Edit": i})
                summary_primers.append({"#Edit": i})
                continue
            all_nicking = []
            all_alternate = []
            for j, p in enumerate(pegRNAs):
                p['#pegRNA'] = j+1
                p['pegRNA'] = p['pegrna']
                p['nsgRNA'] = ''
                try:
                    p['offtargets'] = ','.join(str(e) for e in p['offtargets'][0])
                except KeyError:
                    p['offtargets'] = 'unknown'
                oligos = p.pop('oligos')
                for k in oligos:
                    if isinstance(oligos[k], dict):
                        for k2 in oligos[k]:
                            p[f'{k}_oligo_{k2}'] = oligos[k][k2]
                    else:
                        p[f'{k}_oligo'] = oligos[k]

                nicking = p.pop('nicking')
                alternate_extensions = p.pop('alternate_extensions')
                pegRNAs[j] = {k: p[k] for k in pegRNA_header}
                if cloning_strategy.can_design_nicking:
                    try:
                        for k, v in zip(cloning_strategy.excel_nicking_headers, nicking[0].keys()):
                            pegRNAs[j][k] = nicking[0][v]
                        pegRNAs[j]['nsgRNA'] = nicking[0]['spacer']
                    except IndexError:
                        for k in cloning_strategy.excel_nicking_headers:
                            pegRNAs[j][k] = ''
                    for n in nicking:
                        n['#pegRNA'] = j + 1
                        try:
                            n['offtargets'] = ','.join(str(e) for e in n['offtargets'][0])
                        except KeyError:
                            n['offtargets'] = 'unknown'
                        for k, v in zip(cloning_strategy.excel_nicking_headers, list(nicking[0].keys())):
                            n[k] = nicking[0][v]

                        n = {k: n[k] for k in ['#pegRNA', 'kind', 'position', 'spacer', 'score', 'wt_score', 'offtargets'] \
                                                 + cloning_strategy.excel_nicking_headers}
                        all_nicking.append(n)

                for ext in alternate_extensions:
                    oligos = ext.pop('oligos')
                    ext['#pegRNA'] = j + 1
                    for k, v in oligos.items():
                        ext[f'{k}_oligo'] = v

                    ext = {k: ext[k] for k in ['#pegRNA', 'pbs_length', 'rt_template_length', 'sequence', 'pbs_gc',
                                               'rt_gc'] +
                           cloning_strategy.excel_extension_headers}
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
                    summary_dict = {'#Edit': i, 'LEFT_SEQUENCE': '', 'RIGHT_SEQUENCE': ''}
                    summary_dict.update({k: d[k] for k in sorted(d, key=lambda x: (x.startswith('p'), x.startswith('r')))})
                    summary_primers.append(summary_dict)

            start_row = 1

            column_widths = defaultdict(lambda: 0)

            ws = wb.add_worksheet(wsname)
            writer.sheets[wsname] = ws
            for j, df in enumerate([pegRNAs, all_nicking, all_primers, all_alternate]):
                ws.write(f'A{start_row}', ['pegRNAs', 'nsgRNAs', 'primers', 'alternate extensions'][j], heading)
                df = pd.DataFrame(df)
                for col_idx, column in enumerate(df):
                    column_length = max(df[column].astype(str).map(len).max(), len(column))
                    column_widths[col_idx] = max(column_length, column_widths[col_idx])
                df.to_excel(writer, wsname, index=False, startrow=start_row)
                start_row += len(df.index) + 3
            summary_dict = {'#Edit': i}
            summary_dict.update(pegRNAs[0])
            del summary_dict["#pegRNA"]
            summary_pegrnas.append(summary_dict)
            for col_idx, column_length in column_widths.items():
                writer.sheets[wsname].set_column(col_idx, col_idx, column_length + 3)
            ws._opt_close()

        ws = wb.get_worksheet_by_name('Summary')
        summary_start_row += len(summary.index) + 3
        for title, table in (('pegRNAs', summary_pegrnas), ('Primers', summary_primers)):
            ws.write(f'A{summary_start_row}', title, heading)
            df = pd.DataFrame(table)
            df.to_excel(writer, 'Summary', index=False, startrow=summary_start_row)
            for col_idx, column in enumerate(df):
                column_length = max(df[column].astype(str).map(len).max(), len(column))
                summary_column_widths[col_idx] = max(column_length, summary_column_widths[col_idx])
            summary_start_row += len(table) + 3

        for col_idx, column_length in summary_column_widths.items():
            writer.sheets['Summary'].set_column(col_idx, col_idx, column_length + 3)
        writer.save()

        self.excel_exported = True

    def design_edit(self, i):
        edit = self.edits[i]
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
            edit, sequence_object, padding = EditSerializer.parse_edit_dict(edit_dict, self.organism.pk,
                                                                            dict(
                                                                                {'nuclease': self.nuclease,
                                                                                 'nuclease_options': self.nuclease_options,
                                                                                 'cloning_strategy': self.cloning_strategy,
                                                                                 'cloning_options': self.cloning_options,
                                                                                 'design_primers': self.design_primers
                                                                                 }, **self.options))
            result.update(edit.create_oligos())

            if len(result['pegRNAs']) < 1:
                result['warning'] = 'No pegRNAs found'

            if not (self.run_bowtie and len(self.edits) <= django.conf.settings.DESIGN_MAX_EDITS_BOWTIE):
                for pegRNA in result['pegRNAs']:
                    pegRNA['offtargets'] = [['-']]
                    for nicking in pegRNA['nicking']:
                        nicking['offtargets'] = [['-']]
                for primers in result['primers']:
                    primers['product_count'] = '-'

            result['start'] -= padding
            result['edit'] = edit.__class__.__name__
            result['sequence_object'] = SequenceObjectSerializer(sequence_object).data
            if sequence_type == 'custom':
                result['sequence'] = ','.join([sequence_object.name, sequence_object.sequence])
        except Exception as e:
            result['warning'] = f'An unknown error occurred, {repr(e)}'

        finally:
            self.save_result(result, i)
            best = None
            if len(result['pegRNAs']) > 0:
                best = json.dumps(result['pegRNAs'][0])
            return {
                'sequence': result['sequence'] if sequence_type != 'custom' else result['sequence'].split(',')[0],
                'sequence_type': sequence_type,
                'edit': edit.__class__.__name__,
                'options': options,
                'pegRNAs': len(result['pegRNAs']),
                'warning': result['warning'],
                'best': best,
            }

    def primer_specificity_check(self):
        all_primers = itertools.chain.from_iterable((result['primers'] for result in self.results))
        all_primers = check_primer_specificity(all_primers, self.organism.assembly, self.jobdir, **self.options)

        prevpair = 0
        count = 0
        result_index = 0
        results = self.results
        result = next(results)
        while len(result['primers']) == 0:
            try:
                result = next(results)
            except StopIteration:
                break
            result_index += 1
        try:
            result['primers'][count]['products'] = set()
        except IndexError:
            pass
        for pair, product in all_primers:
            if prevpair != pair:
                result['primers'][count]['product_count'] = len(result['primers'][count]['products'])
                count += 1
                try:
                    result['primers'][count]['products'] = set()
                except IndexError:
                    try:
                        result['primers'] = sorted(result['primers'], key=lambda x: x['product_count'])
                        self.save_result(result, result_index)
                        result = next(results)
                        result_index += 1
                        while len(result['primers']) == 0:
                            result = next(results)
                            result_index += 1
                        count = 0
                        result['primers'][count]['products'] = set()
                    except StopIteration:
                        break

            result['primers'][count]['products'].add(product)
        self.save_result(result, result_index)

    def sgrna_specificity_check(self):
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
                                                                   self.jobdir)
            for i, result in enumerate(self.results):
                for j, hits in enumerate(zip(counts, binders)):
                    for pegRNA in result['pegRNAs']:
                        if pegRNA['nuclease'] == nuclease and pegRNA['spacer'] == spacers[nuclease][j]:
                            pegRNA['offtargets'] = hits
                        for nicking in pegRNA['nicking']:
                            if pegRNA['nuclease'] == nuclease and nicking['spacer'] == spacers[nuclease][j]:
                                nicking['offtargets'] = hits
                self.save_result(result, i)

    def create_oligos(self):
        """Design pegRNAs for job edits. Saves to jobdir folder"""
        #try:
        self.status = JobStatus.FINDING_PEGRNAS
        self.warning = None
        self.summary = []
        self.save()
        for i in range(len(self.edits)):
            self.summary.append(self.design_edit(i))
            self.save()

        nuclease = NUCLEASES[self.nuclease]
        cloning_strategy = nuclease.cloning_strategies[self.cloning_strategy]
        for i, res in enumerate(cloning_strategy.post_process(self.results)):
            self.save_result(res, i)

        if self.run_bowtie:
            if self.design_primers:
                self.status = JobStatus.CHECKING_PRIMER_SPECIFICITY
                self.save()
                self.primer_specificity_check()

            self.status = JobStatus.CHECKING_SGRNA_SPECIFICITY
            self.save()
            self.sgrna_specificity_check()

        self.status = JobStatus.EXPORTING_EXCEL
        self.save()
        self.export_excel()
        self.status = JobStatus.COMPLETED_STATUS
        self.save()
        #except Exception as e:
        #    self.status = JobStatus.FAILED_STATUS
        #    self.warning = f'An unknown error occurred, {repr(e)}'
        #    self.save()

        #    return e


class JobSerializer(serializers.Serializer):
    organism = OrganismSerializer()
    status = serializers.CharField()
    job_name = serializers.CharField()
    job_id = serializers.CharField()
    summary = serializers.ListField(child=serializers.DictField(child=serializers.CharField()))
    edits = serializers.ListField(child=serializers.DictField(child=serializers.CharField()))
    options = serializers.DictField(child=serializers.IntegerField())
    nuclease = serializers.CharField(required=False)
    nuclease_options = serializers.DictField(child=serializers.CharField())
    cloning_strategy = serializers.CharField(required=False)
    cloning_options = serializers.DictField(child=serializers.CharField())
    warning = serializers.CharField()
    run_bowtie = serializers.BooleanField()
    design_primers = serializers.BooleanField()
    excel_exported = serializers.BooleanField()


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


@celery.shared_task()
def create_oligos_background(job_id: str) -> None:
    """Used by celery to run create_oligos in background."""
    job = Job.load_from_disk(job_id)
    e = job.create_oligos()

    if isinstance(e, Exception):
        job.status = JobStatus.FAILED_STATUS
        job.save()
        raise e


@celery.shared_task()
def design_edit_background(edit_index: list, job_id: str):
    #print(indices, job_id)
    job = Job.load_from_disk(job_id)
    results = []
    #for edit_index in indices:
    results.append(job.design_edit(edit_index))
    return results


@celery.shared_task()
def update_summary(results: list, job_id: str):
    if type(results[0]) == list:
        results = list(itertools.chain.from_iterable(results))
    job = Job.load_from_disk(job_id)
    job.summary.extend(results)
    job.save()


@celery.shared_task()
def post_process_job_background(job_id: str):
    job = Job.load_from_disk(job_id)
    nuclease = NUCLEASES[job.nuclease]
    cloning_strategy = nuclease.cloning_strategies[job.cloning_strategy]
    for i, res in enumerate(cloning_strategy.post_process(job.results)):
        job.save_result(res, i)


@celery.shared_task()
def queue_primer_specificity(job_id: str):
    job = Job.load_from_disk(job_id)
    job.status = JobStatus.QUEUED_PRIMER_SPECIFITY
    job.save()


@celery.shared_task()
def primer_specificity_background(job_id: str):
    job = Job.load_from_disk(job_id)
    job.status = JobStatus.CHECKING_PRIMER_SPECIFICITY
    job.save()
    job.primer_specificity_check()
    job.save()

    return job_id


@celery.shared_task()
def queue_spacer_specificity(job_id: str):
    job = Job.load_from_disk(job_id)
    job.status = JobStatus.QUEUED_SGRNA_SPECIFITY
    job.save()


@celery.shared_task()
def spacer_specificity_background(job_id: str):
    job = Job.load_from_disk(job_id)
    job.status = JobStatus.CHECKING_SGRNA_SPECIFICITY
    job.save()
    job.sgrna_specificity_check()
    job.save()

    return job_id


@celery.shared_task()
def init_job(job_id: str):
    job = Job.load_from_disk(job_id)
    job.warning = None
    job.summary = []
    job.status = JobStatus.FINDING_PEGRNAS
    job.save()


@celery.shared_task()
def export_excel(job_id: str):
    job = Job.load_from_disk(job_id)
    job.status = JobStatus.EXPORTING_EXCEL
    job.save()
    job.export_excel()
    job.status = JobStatus.COMPLETED_STATUS
    job.save()


def create_oligos_chain(job_id: str):
    job = Job.load_from_disk(job_id)
    job.status = JobStatus.QUEUED
    job.save()

    edit_dict = job.edits[0].copy()
    edit_dict.update({'organism': job.organism.id})

    edit, sequence_object, padding = EditSerializer.parse_edit_dict(edit_dict, job.organism.id,
                                                                    dict(
                                                                        {'nuclease': job.nuclease,
                                                                         'nuclease_options': job.nuclease_options,
                                                                         'cloning_strategy': job.cloning_strategy,
                                                                         'cloning_options': job.cloning_options,
                                                                         'design_primers': job.design_primers
                                                                         }, **job.options))
    edit_dict['sequence_object'] = SequenceObjectSerializer(sequence_object).data
    job.save_result(edit_dict, 0)

    indices = list(range(len(job.edits)))
    n = DESIGN_EDIT_GROUP_SIZE
    indices = [indices[i * n:(i + 1) * n] for i in range((len(indices) + n - 1) // n )]

    results = init_job.si(job_id)
    for idx in indices:
        results = results | celery.group([design_edit_background.si(i, job_id) for i in idx]) | update_summary.s(job_id)
    results = results | post_process_job_background.si(job_id)
    if job.run_bowtie and len(job.edits) <= django.conf.settings.DESIGN_MAX_EDITS_BOWTIE:
        if job.design_primers and NUCLEASES[job.nuclease].cloning_strategies[job.cloning_strategy].can_design_primers:
            results = results | queue_primer_specificity.si(job_id)
            results = results | primer_specificity_background.si(job_id)
        results = results | queue_spacer_specificity.si(job_id)
        results = results | spacer_specificity_background.si(job_id)
    results = results | export_excel.si(job_id)
    return results.delay()
