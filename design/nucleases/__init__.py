import abc
import dataclasses
import typing

from django.utils.functional import classproperty

NUCLEASES = {}


class OligoDict(typing.TypedDict):
    top: str
    bottom: str


@dataclasses.dataclass()
class Nuclease(abc.ABC):
    """Abstract Edit class.
    Can be subclassed to make functional nucleases
    """

    _target_motif: str = None  # type: ignore
    _pam_motif: str = None  # type: ignore

    _upstream_length: int = 0
    _downstream_length: int = 0
    _spacer_length: int = None  # type: ignore
    _cut_site_position: int = None  # type: ignore

    _scaffolds: typing.Dict = dataclasses.field(default_factory=dict)
    _default_scaffold: str = None  # type: ignore

    _cloning_strategies: typing.Dict = dataclasses.field(default_factory=dict)
    _design_strategies: typing.Dict = dataclasses.field(default_factory=dict)

    _options: typing.Dict = dataclasses.field(default_factory=dict)

    @classproperty
    def target_motif(cls):
        return cls._target_motif

    @classproperty
    def pam_motif(cls):
        return cls._pam_motif

    @classproperty
    def upstream_length(cls):
        return cls._upstream_length

    @classproperty
    def downstream_length(cls):
        return cls._downstream_length

    @classproperty
    def spacer_length(cls):
        return cls._spacer_length

    @classproperty
    def cut_site_position(cls):
        return cls.upstream_length + cls._cut_site_position

    @classproperty
    def scaffolds(cls):
        try:
            return cls._scaffolds
        except AttributeError:
            return {}

    @classproperty
    def default_scaffold(cls):
        return cls._default_scaffold

    @classproperty
    def cloning_strategies(cls):
        try:
            return cls._cloning_strategies
        except AttributeError:
            return {}

    @classproperty
    def design_strategies(cls):
        try:
            return cls._design_strategies
        except AttributeError:
            return {}

    @classproperty
    def downstream_from_cut_site(cls):
        return (cls.spacer_length + len(cls.pam_motif) + cls.downstream_length) - cls._cut_site_position

    @classmethod
    def scaffold_name_to_sequence(cls, scaffold_name: typing.Optional[str] = None) -> str:
        """Get scaffold from name"""
        if scaffold_name is None:
            scaffold = cls.scaffolds[cls.default_scaffold]
        else:
            scaffold = cls.scaffolds[scaffold_name]

        return scaffold

    @classmethod
    def get_design_strategy(cls, strategy):
        if strategy is None:
            strategy = next(iter(cls.design_strategies.keys()))
        return cls.design_strategies[strategy]

    @classmethod
    def score_spacers(cls, spacers: typing.List) -> typing.List:
        """Method to calculate score for a list of spacers."""
        return [0] * len(spacers)

    @classmethod
    @abc.abstractmethod
    def _filter_offtarget(cls, seq):
        """Should return True, if seq is a valid off target"""

    @classmethod
    def filter_extension(cls, seq, **options):
        return False

    @classmethod
    def get_cloning_strategy(cls, strategy):
        if isinstance(strategy, str):
            strategy = cls.cloning_strategies[strategy]
        return strategy

    @classproperty
    def options(cls):
        try:
            return cls._options
        except AttributeError:
            return dict()


@dataclasses.dataclass()
class BaseCloningStrategy(abc.ABC):

    can_design_primers: bool = True
    can_design_nicking: bool = True

    _options: typing.Dict = dataclasses.field(default_factory=dict)

    @classmethod
    def _spacer_to_cloning(cls, spacer_sequence: str) -> str:
        # Overwrite in subclasses to modify spacer sequence for cloning
        # Eg. prepend a 'G' for U6 promoters
        return spacer_sequence

    @classmethod
    @abc.abstractmethod
    def design_cloning(cls, **options):
        pass

    @classmethod
    def post_process(cls, designs):
        return designs

    @classmethod
    def can_express(cls, sequence):
        return True

    @classproperty
    def options(cls):
        try:
            return cls._options
        except AttributeError:
            return dict()


@dataclasses.dataclass()
class BaseDesignStrategy(abc.ABC):

    _options: typing.Dict = dataclasses.field(default_factory=dict)

    @classmethod
    @abc.abstractmethod
    def find_spacers(cls, reference_sequence: str, altered_sequence: str, start: int, end: int,
                     spacer_search_range: int, **options) -> typing.List:
        """
        Method to find pegRNA spacer candidates.
        Should return spacers within range of alteration
        """

    @classmethod
    @abc.abstractmethod
    def find_nicking_spacers(cls, reference_sequence: str, altered_sequence: str, spacer_strand: typing.Literal[1, -1],
                             cut_site: int, scaffold: str, nicking_range: int, cloning_method,
                             **options) -> typing.List:
        """
        Method to find nicking spacers for a given pegRNA (cut site and strand).
        Should return spacers on opposite strand that cut within range
        """

    @staticmethod
    def make_extension_sequence(nuclease, reference_sequence, altered_sequence, spacer_strand, spacer_cut_site,
                                cut_dist,
                                alteration_length, pbs_min_length, pbs_max_length, rt_min_length, rt_max_length,
                                strategy,
                                **options):
        """Method to create the pegRNA extension sequence.
        pegRNA extension sequences consist of a PBS and a RT template.
        The PBS is upstream of the cut site. The RT template is downstream of the cut site and contains the edit sequence.
        """

    @classproperty
    def options(cls):
        try:
            return cls._options
        except AttributeError:
            return dict()
