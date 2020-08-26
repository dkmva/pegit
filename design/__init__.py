from django.conf import settings
from django.core.exceptions import ImproperlyConfigured

from design.edits import load_edits
from design.nucleases import load_nucleases

load_edits()
load_nucleases()


try:
    DESIGN_TWO_BIT_TO_FA_PATH = settings.DESIGN_TWO_BIT_TO_FA_PATH
except AttributeError:
    raise ImproperlyConfigured('DESIGN_TWO_BIT_TO_FA_PATH must be set')

try:
    DESIGN_PRIMER3_PATH = settings.DESIGN_PRIMER3_PATH
except AttributeError:
    raise ImproperlyConfigured('DESIGN_PRIMER3_PATH must be set')

try:
    DESIGN_PRIMER3_CONFIG_PATH = settings.DESIGN_PRIMER3_CONFIG_PATH
except AttributeError:
    raise ImproperlyConfigured('DESIGN_PRIMER3_CONFIG_PATH must be set')

try:
    DESIGN_BOWTIE_PATH = settings.DESIGN_BOWTIE_PATH
except AttributeError:
    raise ImproperlyConfigured('DESIGN_BOWTIE_PATH must be set')
try:
    DESIGN_ASSEMBLIES_FOLDER = settings.DESIGN_ASSEMBLIES_FOLDER
except AttributeError:
    raise ImproperlyConfigured('DESIGN_ASSEMBLIES_FOLDER must be set')

try:
    DESIGN_BOWTIE_GENOMES_FOLDER = settings.DESIGN_BOWTIE_GENOMES_FOLDER
except AttributeError:
    raise ImproperlyConfigured('DESIGN_BOWTIE_GENOMES_FOLDER must be set')

try:
    DESIGN_CODON_USAGE_FOLDER = settings.DESIGN_CODON_USAGE_FOLDER
except AttributeError:
    raise ImproperlyConfigured('DESIGN_CODON_USAGE_FOLDER must be set')

try:
    DESIGN_OUTPUT_FOLDER = settings.DESIGN_OUTPUT_FOLDER
except AttributeError:
    raise ImproperlyConfigured('DESIGN_OUTPUT_FOLDER must be set')

try:
    DESIGN_CLINVAR_ORGANISM = settings.DESIGN_CLINVAR_ORGANISM
except AttributeError:
    raise ImproperlyConfigured('DESIGN_CLINVAR_ORGANISM must be set')

try:
    DESIGN_ENTREZ_EMAIL = settings.DESIGN_ENTREZ_EMAIL
except AttributeError:
    raise ImproperlyConfigured('DESIGN_ENTREZ_EMAIL must be set')

try:
    DESIGN_BOWTIE_THREADS = settings.DESIGN_BOWTIE_THREADS
except AttributeError:
    settings.DESIGN_BOWTIE_THREADS = 1
