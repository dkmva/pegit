import argparse
import json
import os
import pathlib
import subprocess

from Bio import SeqIO
import yaml

import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'prime.settings')
django.setup()

from design.interface import Job
from design.models import Organism


def list_organisms(namespace) -> None:
    print('Available Organisms')
    print('-------------------')
    for organism in Organism.objects.all():
        print(organism.pk, organism)


def add_organism(namespace) -> None:
    options = vars(namespace)
    name = options['name']
    assembly = options['assembly']
    gff = options['gff']
    source = options['source']
    kind = options['kind']
    codon_table = options['codon_table']
    sequence_search = options['sequence_search']
    annotation_name = options['annotation_name']
    scaffolds = options['scaffolds']
    if scaffolds is None:
        scaffolds = ''
    else:
        with open(scaffolds) as f:
            scaffolds = f.read()

    Organism.add_to_database(name=name, assembly=assembly, file_path=gff, source=source, kind=kind,
                             codon_table=codon_table, sequence_search=sequence_search, annotation_name=annotation_name,
                             scaffolds=scaffolds)


def run_organism(namespace) -> None:
    options = vars(namespace)
    del options['func']
    organism = options.pop('organism')
    edits = options.pop('edits')
    output_folder = options.pop('output_folder')
    job_name = options.pop('job_name')
    nuclease = options.pop('nuclease')
    nuclease_options = options.pop('nuclease_options')
    cloning_strategy = options.pop('cloning_strategy')
    cloning_options = options.pop('cloning_options')
    run_bowtie = options.pop('run_bowtie')
    design_primers = options.pop('design_primers')
    organism = Organism.objects.get(pk=organism)

    j = Job(organism, options=options, job_id=job_name, job_name=job_name, output_folder=output_folder,
            nuclease=nuclease, nuclease_options=nuclease_options,
            cloning_strategy=cloning_strategy, cloning_options=cloning_options,
            run_bowtie=run_bowtie, design_primers=design_primers)
    print('Parsing edit list')
    j.import_edit_list(edits)
    print('Designing oligos')
    j.create_oligos()
    print(f'Oligos saved to {j.jobdir}')


def run_clinvar(namespace) -> None:
    options = vars(namespace)
    del options['func']
    organism = Organism.objects.get(assembly=django.conf.settings.DESIGN_CLINVAR_ORGANISM)
    edits = options.pop('edits')
    output_folder = options.pop('output_folder')
    job_name = options.pop('job_name')
    nuclease = options.pop('nuclease')
    nuclease_options = options.pop('nuclease_options')
    cloning_strategy = options.pop('cloning_strategy')
    cloning_options = options.pop('cloning_options')
    run_bowtie = options.pop('run_bowtie')
    design_primers = options.pop('design_primers')

    j = Job(organism, options=options, job_id=job_name, job_name=job_name, output_folder=output_folder,
            nuclease=nuclease, nuclease_options=nuclease_options,
            cloning_strategy=cloning_strategy, cloning_options=cloning_options,
            run_bowtie=run_bowtie, design_primers=design_primers)
    print('Parsing edit list')
    j.import_clinvar_list(edits)
    print('Designing oligos')
    j.create_oligos()
    print(f'Oligos saved to {j.jobdir}')


def _make_2bit_and_scaffold_list(fasta, name) -> str:
    if name is None:
        name = pathlib.Path(fasta).stem

    print('Running faToTwoBit')
    subprocess.run(
        f'{os.path.join(os.path.dirname(django.conf.settings.DESIGN_TWO_BIT_TO_FA_PATH), "faToTwoBit")} "{fasta}" "{django.conf.settings.DESIGN_ASSEMBLIES_FOLDER}/{name}.2bit"',
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='ascii', check=True)
    

    print('Getting record ids')
    scaffolds = [record.id for record in SeqIO.parse(fasta, 'fasta')]
    scaffolds_path = f'{django.conf.settings.DESIGN_ASSEMBLIES_FOLDER}/{name}'
    with open(scaffolds_path, 'w') as f:
        f.write(' '.join(scaffolds))
    print('done')

    return scaffolds_path


def make_2bit_and_scaffold_list(namespace) -> None:
    options = vars(namespace)
    fasta = options['fasta']
    name = options['name']
    _make_2bit_and_scaffold_list(fasta, name)


def _make_bowtie(fasta, name) -> None:
    if name is None:
        name = pathlib.Path(fasta).stem

    print('Building bowtie index')
    try:
        os.makedirs(os.path.join(django.conf.settings.DESIGN_BOWTIE_GENOMES_FOLDER, name))
    except FileExistsError:
        pass
    subprocess.run(
        f'{django.conf.settings.DESIGN_BOWTIE_PATH}-build "{fasta}" "{django.conf.settings.DESIGN_BOWTIE_GENOMES_FOLDER}/{name}/{name}"',
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='ascii', check=True)

    print('done')


def make_bowtie(namespace) -> None:
    options = vars(namespace)
    fasta = options['fasta']
    name = options['name']

    _make_bowtie(fasta, name)


def make_2bit_and_bowtie(namespace) -> None:
    make_2bit_and_scaffold_list(namespace)
    make_bowtie(namespace)


def add_organism_folder(namespace) -> None:
    options = vars(namespace)
    path = options['path']
    with open(os.path.join(path, 'info.yaml')) as f:
        conf = yaml.safe_load(f)

    assembly = conf['assembly']
    fasta = os.path.join(path, 'assembly.fa')
    gff = os.path.join(path, 'genes.gff3')
    scaffolds = _make_2bit_and_scaffold_list(fasta, assembly)
    _make_bowtie(fasta, assembly)

    with open(scaffolds) as f:
        scaffolds = f.read()

    Organism.add_to_database(**conf, file_path=gff, scaffolds=scaffolds)



def main():
    parser = argparse.ArgumentParser('pegIT')
    parser.set_defaults(func=lambda x: parser.print_help())
    subparsers = parser.add_subparsers()

    # Functions related to Organisms model
    organisms_parser = subparsers.add_parser('organisms')
    organisms_parser.set_defaults(func=lambda x: organisms_parser.print_usage())
    organisms_subparsers = organisms_parser.add_subparsers()

    # Print organisms in the database
    list_parser = organisms_subparsers.add_parser('list')
    list_parser.set_defaults(func=list_organisms)

    # Add new organism to the database
    add_parser = organisms_subparsers.add_parser('add')
    add_parser.set_defaults(func=add_organism)
    add_parser.add_argument('--name', required=True, help='Name of the organism')
    add_parser.add_argument('--assembly', required=True, help='Name of the 2bit genome file')
    add_parser.add_argument('--gff', required=True, help='GFF3 file with annotations')
    add_parser.add_argument('--kind', required=True, choices=['gencode', 'refseq', 'ensembl'], help='Kind of GFF3 file')
    add_parser.add_argument('--source', required=True, help='URL to the genome')
    add_parser.add_argument('--codon_table', required=False, help='Name of codon usage json file, defaults to human codon usage if not specified')
    add_parser.add_argument('--sequence_search', required=True, help='URL to sequence search')
    add_parser.add_argument('--annotation_name', required=False, help='Name of the annotation, defaults to the basename of the GFF3 file.')
    add_parser.add_argument('--scaffolds', required=False, help='path to file with scaffolds list')

    # Add new organism to the database
    add_parser = organisms_subparsers.add_parser('folder')
    add_parser.set_defaults(func=add_organism_folder)
    add_parser.add_argument('path', help='path to folder')

    # Functions for designing pegRNAs
    design_parser = subparsers.add_parser('design')
    design_parser.set_defaults(func=lambda x: design_parser.print_usage())
    design_subparsers = design_parser.add_subparsers()

    # Organism mode
    organism_parser = design_subparsers.add_parser('organism')
    organism_parser.set_defaults(func=run_organism)
    organism_parser.add_argument('organism')
    organism_parser.add_argument('edits', help='Edit list')
    organism_parser.add_argument('--output_folder', help='Output folder', default=os.getcwd())
    organism_parser.add_argument('--job_name', help='Project name', default='pegIT')
    organism_parser.add_argument('--nuclease', help='Nuclease', default='SpCas9')
    organism_parser.add_argument('--cloning_strategy', help='Cloning plasmid', default='pegRNA-GG-acceptor')
    organism_parser.add_argument('--nuclease_options', help='dictionary with nuclease options', type=json.loads)
    organism_parser.add_argument('--cloning_options', help='dictionary with cloning options', type=json.loads)
    organism_parser.add_argument('--run_bowtie', dest='run_bowtie', action='store_true')
    organism_parser.add_argument('--no_run_bowtie', dest='run_bowtie', action='store_false')
    organism_parser.set_defaults(run_bowtie=True)
    organism_parser.add_argument('--design_primers', dest='design_primers', action='store_true')
    organism_parser.add_argument('--no_design_primers', dest='design_primers', action='store_false')
    organism_parser.set_defaults(design_primers=True)
    for opt, val in django.conf.settings.DESIGN_CONF['default_options'].items():
        organism_parser.add_argument(f'-{opt}', required=False, default=val)

    # ClinVar mode
    clinvar_parser = design_subparsers.add_parser('clinvar')
    clinvar_parser.set_defaults(func=run_clinvar)
    clinvar_parser.add_argument('edits', help='Edit list')
    clinvar_parser.add_argument('--output_folder', help='Output folder', default=os.getcwd())
    clinvar_parser.add_argument('--job_name', help='Project name', default='pegIT')
    clinvar_parser.add_argument('--nuclease', help='Nuclease', default='SpCas9')
    clinvar_parser.add_argument('--cloning_strategy', help='Cloning plasmid', default='pegRNA-GG-acceptor')
    clinvar_parser.add_argument('--nuclease_options', help='dictionary with nuclease options', type=json.loads)
    clinvar_parser.add_argument('--cloning_options', help='dictionary with cloning options', type=json.loads)
    clinvar_parser.add_argument('--run_bowtie', dest='run_bowtie', action='store_true')
    clinvar_parser.add_argument('--no_run_bowtie', dest='run_bowtie', action='store_false')
    clinvar_parser.set_defaults(run_bowtie=True)
    clinvar_parser.add_argument('--design_primers', dest='design_primers', action='store_true')
    clinvar_parser.add_argument('--no_design_primers', dest='design_primers', action='store_false')
    clinvar_parser.set_defaults(design_primers=True)
    for opt, val in django.conf.settings.DESIGN_CONF['default_options'].items():
        clinvar_parser.add_argument(f'-{opt}', required=False, default=val)


    # utils
    utils_parser = subparsers.add_parser('utils')
    utils_parser.set_defaults(func=lambda x: utils_parser.print_usage())
    utils_subparsers = utils_parser.add_subparsers()

    # Make 2bit file and get scaffolds from fasta
    mk2bit_parser = utils_subparsers.add_parser('2bit')
    mk2bit_parser.set_defaults(func=make_2bit_and_scaffold_list)

    mk2bit_parser.add_argument('fasta', help='fasta file')
    mk2bit_parser.add_argument('--name', help='assembly name', required=False)

    # Build bowtie index
    bowtie_parser = utils_subparsers.add_parser('bowtie')
    bowtie_parser.set_defaults(func=make_bowtie)

    bowtie_parser.add_argument('fasta', help='fasta file')
    bowtie_parser.add_argument('--name', help='assembly name', required=False)

    # Build 2bit and bowtie index
    bowtie_parser = utils_subparsers.add_parser('both')
    bowtie_parser.set_defaults(func=make_2bit_and_bowtie)

    bowtie_parser.add_argument('fasta', help='fasta file')
    bowtie_parser.add_argument('--name', help='assembly name', required=False)

    args = parser.parse_args()
    args.func(args)
    return


if __name__ == '__main__':
    main()
