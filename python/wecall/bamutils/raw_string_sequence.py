# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.read_sequence import ReadSequence, ReadSequenceWithCoverage
from wecall.bamutils.sequence import Sequence
from wecall.bamutils.sequence_position import DELETED_BASE
from wecall.bamutils.sequence_quality import SequenceQuality
from wecall.common.exceptions import EchidnaException
import re

from wecall.genomics.reference_chromosome import ReferenceChromosome


class RawStringSequences(object):

    def __init__(self):
        self.reference_string = ""
        self.__sequence_string = ""
        self.quality_string = ""
        self.pos_from = None
        self.is_ongoing = False

    @property
    def sequence_string(self):
        if not self.is_forward_seq() and not self.is_reverse_seq():
            raise EchidnaException(
                "Illegal character in sequence {!r}".format(
                    self.__sequence_string))

        return self.__sequence_string

    def add_position(self, sequence_position, ref_pos):
        self.reference_string += sequence_position.ref_char
        self.__sequence_string += sequence_position.seq_char

        if sequence_position.seq_char != DELETED_BASE:
            self.quality_string += sequence_position.qual_char

        self.is_ongoing = True

        if self.pos_from is None:
            self.pos_from = ref_pos

    def is_forward_seq(self):
        return re.match(r'^[ACGTURYKMSWBDHVN\*\.]*\Z', self.__sequence_string)

    def is_reverse_seq(self):
        return re.match(r'^[acgturykmswbdhvn\*\,]*\Z', self.__sequence_string)

    def build_annotated_seq(
            self, n_fwd, n_rev, mapping_quality, insert_size,
            read_id, read_flags, cigar_string, read_start, read_mate_start
    ):
        reference = ReferenceChromosome(self.reference_string, self.pos_from)
        sequence = Sequence(
            reference, self.sequence_string.replace(
                ",", ".").upper(), cigar_string)
        quality = SequenceQuality(self.quality_string)

        read_sequence = ReadSequence(
            sequence,
            quality,
            mapping_quality,
            insert_size,
            read_id,
            read_flags,
            read_start,
            read_mate_start)
        if n_fwd is not None:
            return [ReadSequenceWithCoverage(read_sequence, n_fwd, n_rev)]
        elif self.is_reverse_seq():
            return [ReadSequenceWithCoverage(read_sequence, 0, 1)]
        elif self.is_forward_seq():
            return [ReadSequenceWithCoverage(read_sequence, 1, 0)]
        else:
            raise EchidnaException(
                "Raw sequence: {} is neither forward or reverse".format(self))
