import abc
import dataclasses


@dataclasses.dataclass()
class BaseCloningStrategy(abc.ABC):

    can_design_primers: bool = True
    can_design_nicking: bool = True
    allow_extension_filtering: bool = True

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
