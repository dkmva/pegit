"""
prime_design edits
"""
import abc
import typing

import django
import inflection

from design.oligos import AlterationTracker

EDITS = {}


class AbstractEdit(abc.ABC):
    """Abstract Edit class.
   Can be subclassed to make functional edits, by overriding .run().
   """

    # Default options
    default_options = {}

    # Required options
    required = ()

    # Set to true, if allowed outside coding sequences.
    allow_genomic = False

    def __init__(self, sequence_object, option_string, **options):
        self.options = {}
        self.sequence_object = sequence_object
        options.update(self.parse_option_string(option_string))
        self.repair = 'repair' in options and options['repair'] == 'true'
        self.silence_pam = options.get('silence_pam', False)
        if self.silence_pam == 'false':
            self.silence_pam = False
        try:
            self.pbs_length = int(options['pbs_length'])
        except (ValueError, KeyError):
            self.pbs_length = 'auto'
        try:
            self.rt_template_length = int(options['rt_template_length'])
        except (ValueError, KeyError):
            self.rt_template_length = 'auto'
        self.update_options(options)

    @abc.abstractmethod
    def run(self) -> AlterationTracker:
        """
        Run method. Must be overridden in subclasses.

        Should call validate options and return an AlterationTracker object
        """

    def get_tracker(self) -> AlterationTracker:
        tracker = AlterationTracker(
            self.sequence_object.upstream + self.sequence_object.sequence + self.sequence_object.downstream)
        offset = len(self.sequence_object.upstream)
        for dgn in [(dgn[0], dgn[1]+offset) for dgn in self.sequence_object.degenerate_sequence]:
            tracker.add_degenerate(*dgn)

        return tracker

    def create_oligos(self):
        tracker = self.run()
        return tracker.make_oligos(repair=self.repair, pbs_length=self.pbs_length, rt_template_length=self.rt_template_length, silence_pam=self.silence_pam, **self.options)

    @staticmethod
    def parse_option_string(option_string):
        options = option_string.split(',')
        parsed_options = {}
        for option in options:
            try:
                k, v = option.split('=')
                parsed_options[inflection.underscore(k)] = v
            except ValueError:
                raise ValueError(f"Invalid option '{option}'")
        return parsed_options

    @classmethod
    def validate_options(cls, options: dict, sequence_object):
        """
        Validate options against a sequence object. Default checks whether required options are supplied.
        Should be overriden in subclasses and called by .run()

        :param options: dictionary with edit options
        :param sequence_object: sequence_object to validate against
        :return: None
        """
        for option in cls.default_options:
            if option in cls.required and option not in options:
                raise ValueError(f'Required option {option} not supplied')

    def update_options(self, options: dict) -> None:
        """Add options to instance. Adds class specific and general default options if not supplied."""
        self.validate_options(options, self.sequence_object)
        for option in self.default_options:
            if option in options:
                self.options[option] = self.default_options[option][0](options[option])
            else:
                self.options[option] = self.default_options[option][1]
        defaults = django.conf.settings.DESIGN_CONF['default_options']
        for option in defaults:
            if option in options:
                self.options[option] = int(options[option])
            else:
                self.options[option] = defaults[option]


def register_edit(edit: typing.Type[AbstractEdit]) -> None:
    """Register an edit with the global EDITS dict."""
    EDITS[edit.__name__] = edit


def load_edits():
    for module in django.conf.settings.DESIGN_CONF['edits']:
        __import__(module, globals(), locals())
