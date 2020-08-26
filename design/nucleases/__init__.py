import abc
import dataclasses
import typing

import django


NUCLEASES = {}


class OligoDict(typing.TypedDict):
    top: str
    bottom: str


@dataclasses.dataclass()
class Nuclease(abc.ABC):
    """Abstract Edit class.
   Can be subclassed to make functional nucleases
    """

    target_motif: str = None  # type: ignore
    pam_motif: str = None  # type: ignore
    cut_site_position: int = None  # type: ignore
    downstream_from_cut_site: int = None  # type: ignore

    scaffolds: typing.Dict = dataclasses.field(default_factory=dict)
    default_scaffold: str = None  # type: ignore

    @classmethod
    def _scaffold_name_to_sequence(cls, scaffold_name: typing.Optional[str] = None) -> str:
        """Get scaffold from name"""
        if scaffold_name is None:
            scaffold = cls.scaffolds[cls.default_scaffold]
        else:
            scaffold = cls.scaffolds[scaffold_name]

        return scaffold

    @classmethod
    def _spacer_to_cloning(cls, spacer_sequence: str, scaffold_name: typing.Optional[str] = None) -> typing.Tuple[str, str]:
        scaffold = cls._scaffold_name_to_sequence(scaffold_name)

        target = spacer_sequence.upper()
        if not target.startswith('G'):
            target = 'g' + target

        return target, scaffold

    @classmethod
    def score_spacers(cls, spacers: typing.List) -> typing.List:
        """Method to calculate score for a list of spacers."""
        return [None] * len(spacers)

    @classmethod
    @abc.abstractmethod
    def make_scaffold_oligos(cls, scaffold_name: str = None) -> OligoDict:
        """
        Method to convert a scaffold name to oligos.
        Use cls._scaffold_name_to_sequence(scaffold_name) to get scaffold sequence from name.
        """

    @classmethod
    @abc.abstractmethod
    def make_spacer_oligos(cls, spacer_sequence: str, scaffold_name: typing.Optional[str] = None) -> OligoDict:
        """
        Method to convert a spacer sequence, and optionally a scaffold name to pegRNA spacer oligos.
        Use cls._scaffold_name_to_sequence(scaffold_name) to get scaffold sequence from name.
        """

    @classmethod
    @abc.abstractmethod
    def make_nicking_oligos(cls, spacer_sequence: str, scaffold_name: typing.Optional[str] = None) -> OligoDict:
        """
        Method to convert a spacer sequence, and optionally a scaffold name to nicking oligos.
        Use cls._scaffold_name_to_sequence(scaffold_name) to get scaffold sequence from name.
        """

    @classmethod
    @abc.abstractmethod
    def find_spacers(cls, reference_sequence: str, mutated_sequence: str, start: int, end: int,
                     spacer_search_range: int, **options) -> typing.List:
        """
        Method to find pegRNA spacer candidates.
        Should return spacers within range of mutation
        """

    @classmethod
    @abc.abstractmethod
    def find_nicking_spacers(cls, reference_sequence: str, mutated_sequence: str, spacer_strand: typing.Literal[1, -1],
                             cut_site: int, scaffold: str, nicking_range: int, **options) -> typing.List:
        """
        Method to find nicking spacers for a given pegRNA (cut site and strand).
        Should return spacers on opposite strand that cut within range
        """

    @classmethod
    @abc.abstractmethod
    def _make_pbs_sequence(cls, reference: str, pbs_min_length: int, pbs_max_length: int,
                           **options) -> typing.Tuple[str, typing.List]:
        """
        Method to find PBS sequence for a pegRNA spacer.
        Should return optimal sequence and a list of alternatives
        """

    @classmethod
    @abc.abstractmethod
    def _make_rt_sequence(cls, reference, cut_dist, nucleotide_difference, mutation_length, rt_min_length, rt_max_length, **options) -> typing.Tuple[str, typing.List]:
        """
        Method to find RT template sequence for a pegRNA spacer.
        Should return optimal sequence and a list of alternatives
        """


def register_nuclease(nuclease):
    """Register a nuclease with the global NUCLEASES dict."""
    NUCLEASES[nuclease.__name__] = nuclease


def load_nucleases():
    for module in django.conf.settings.DESIGN_CONF['nucleases']:
        __import__(module, globals(), locals())

