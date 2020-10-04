import importlib
import inspect

from django.conf import settings
from django.core.exceptions import ImproperlyConfigured

from design.edits import AbstractEdit, EDITS
from design.nucleases import Nuclease, NUCLEASES


def auto_import(module_list, base_class, container):
    """Dynamic import of nucleases and edits."""
    for module in module_list:
        module = importlib.import_module(module)
        for attribute in dir(module):
            attribute = getattr(module, attribute)
            if inspect.isclass(attribute) and issubclass(attribute, base_class):
                # Add all fully implemented to dict
                if len(attribute.__dict__['__abstractmethods__']) == 0:
                    container[attribute.__name__] = attribute


auto_import(settings.DESIGN_CONF['edits'], AbstractEdit, EDITS)
auto_import(settings.DESIGN_CONF['nucleases'], Nuclease, NUCLEASES)
NUCLEASES = dict({settings.DESIGN_CONF['default_nuclease']: NUCLEASES[settings.DESIGN_CONF['default_nuclease']]},  **{k:v for (k,v) in NUCLEASES.items()})


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
