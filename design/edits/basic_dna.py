"""
Module for basic DNA alterations
"""
import itertools

from . import AbstractEdit

from design.helpers import is_valid_dna


class Insertion(AbstractEdit):
    """Insert nucleotides.

    ###### Options
    * *insert (str)*: String of nucleotides to insert.
    * *position (int)*: Position of insertion.

    ###### Usage
    DNA sequence is: ATGCGAATATTG
    * insert=ATA, position=7
        * New DNA sequence is: ATGCGAataATATTG
    """
    default_options = {'insert': (str, '', '[ACGTacgt]+'),
                       'position': (int, None, None)}
    required = ('insert', 'position')
    allow_genomic = True

    def run(self):
        insert, position = self.validate_options(self.options, self.sequence_object)
        position = self.sequence_object.contextualize_alteration(position)
        tracker = self.get_tracker()
        for i, nt in enumerate(insert):
            tracker.insert(nt, position+i, force=True)
        return tracker

    @classmethod
    def validate_options(cls, options, sequence_object):
        insert = options['insert']
        position = int(options['position']) - 1
        super().validate_options(options, sequence_object)
        if not is_valid_dna(insert):
            raise ValueError(f'Insert ({insert}) is not valid DNA.')

        return insert, position


class Deletion(AbstractEdit):
    """Delete nucleotides.

    ###### Options
    * *delete (str / int)*: What to delete. Can be either a string of nucleotides, or number of nucleotides.
    * *position (int)*: Position to delete from.

    ###### Usage
    DNA sequence is: ATGCGAATATTG
    * delete=ATA, position=7
        * New DNA sequence is: ATGCGATTG

    * delete=6, position=4
        * New DNA sequence is: ATGTTG
    """
    default_options = {'delete': (str, '', '[ACGTacgt]+|\d+'),
                       'position': (int, None, None)}
    required = ('delete', 'position')
    allow_genomic = True

    def run(self):
        to_delete, position = self.validate_options(self.options, self.sequence_object)
        position = self.sequence_object.contextualize_alteration(position)
        tracker = self.get_tracker()
        for _ in enumerate(to_delete):
            tracker.delete(position, force=True)
        return tracker

    @classmethod
    def validate_options(cls, options, sequence_object):
        super().validate_options(options, sequence_object)
        deletion = options['delete']
        position = int(options['position']) - 1
        try:
            deletion_length = int(deletion)
            to_delete = sequence_object.sequence[position:position + deletion_length]
        except ValueError:
            to_delete = deletion
            deletion_length = len(to_delete)

        target = sequence_object.sequence[position:position + deletion_length]
        if target.upper() != to_delete.upper():
            raise ValueError(f"Tried to delete sequence '{to_delete}', but '{target}', was found in the sequence")

        return to_delete, position


class Substitution(AbstractEdit):
    """Substitute nucleotides.

    ###### Options
    * *from (str)*: Nucleotides to replace.
    * *to (str)*: Nucleotides to replace with.
    * *position (int)*: Position of substitution.

    ###### Usage
    DNA sequence is: ATGCGAATATTG
    * from=ATA, to=GGC, position=7
        * New DNA sequence is: ATGCGAggcTTG

    * from=ATA, to=GGCACG, position=7
        * New DNA sequence is: ATGCGAggcacgTTG

    """
    default_options = {'from': (str, '', '[ACGTacgt]+'),
                       'to': (str, '', '[ACGTacgt]+'),
                       'position': (int, None, None)}
    required = ('from', 'to', 'position')
    allow_genomic = True

    def run(self):
        from_, to_, position = self.validate_options(self.options, self.sequence_object)
        position = self.sequence_object.contextualize_alteration(position)
        tracker = self.get_tracker()

        deletions = 0
        for i, (wt, alt) in enumerate(itertools.zip_longest(from_, to_, fillvalue='-'),
                                      position):
            if wt == '-':
                tracker.insert(alt, i - deletions, force=True)
            elif alt == '-':
                assert tracker[i - deletions].lower() == wt.lower()
                tracker.delete(i - deletions, force=True)
                deletions += 1
            else:
                assert tracker[i - deletions].lower() == wt.lower()
                tracker.substitute(alt, i - deletions, force=True)

        return tracker

    @classmethod
    def validate_options(cls, options, sequence_object):
        super().validate_options(options, sequence_object)
        from_ = options['from']
        to_ = options['to']
        if not is_valid_dna(from_):
            raise ValueError(f'From ({from_}) is not valid DNA.')
        if not is_valid_dna(to_):
            raise ValueError(f'To ({to_}) is not valid DNA.')
        position = int(options['position']) - 1
        found = sequence_object.sequence[position:position + len(from_)]
        if found.upper() != from_.upper():
            raise ValueError(f"Tried to substitute sequence '{from_}', but '{found}', was found in the sequence at position {position+1}")
        if from_.upper() == to_.upper():
            raise ValueError(f"Substituting '{to_}' with '{from_}' will not result in any edits")
        return from_, to_, position
