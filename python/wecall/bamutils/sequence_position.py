# All content Copyright (C) 2018 Genomics plc
from wecall.common.exceptions import EchidnaException

MATCHING_BASE = "."
DELETED_BASE = "*"
MISSING_BASE = " "


class SequencePosition(object):

    def __init__(self, ref_char, seq_char, qual_char):
        self.__validate_input(ref_char, seq_char, qual_char)

        self.ref_char = ref_char
        self.seq_char = seq_char
        self.is_gap = self.__calculate_is_gap()
        self.qual_char = qual_char

        self.__validate_character_combination()

    def update_ref_pos(self, ref_pos):
        if self.ref_char != DELETED_BASE:
            return ref_pos + 1
        else:
            return ref_pos

    def __calculate_is_gap(self):
        return self.seq_char == MISSING_BASE

    @staticmethod
    def __validate_input(ref_char, seq_char, qual_char):
        if not all(len(c) == 1 for c in [ref_char, seq_char, qual_char]):
            raise EchidnaException(
                "All characters at sequence position has to be of length 1.")

        if ref_char == MISSING_BASE:
            raise EchidnaException("Missing reference character.")

    def __validate_character_combination(self):
        if self.ref_char == DELETED_BASE and self.seq_char == MATCHING_BASE:
            raise EchidnaException(
                "Invalid character combination: ref char = {}, sequence char = {}".format(
                    self.ref_char, self.seq_char))

        if self.seq_char == DELETED_BASE and self.qual_char != MISSING_BASE:
            raise EchidnaException(
                "Cannot assign base quality to a deleted base.")
        if self.is_gap and self.qual_char != MISSING_BASE:
            raise EchidnaException("Cannot assign base quality inside a gap.")
