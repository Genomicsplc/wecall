# All content Copyright (C) 2018 Genomics plc
import itertools


class Cigar(object):
    MATCH = object()
    INSERTION = object()
    DELETION = object()
    SKIP = object()
    SOFT_CLIP = object()
    HARD_CLIP = object()
    PADDING = object()
    SEQUENCE_MATCH = object()
    SEQUENCE_MISMATCH = object()

    def __init__(self, cigar_sequence=None):
        if cigar_sequence is None:
            cigar_sequence = []
        # sequence of cigar components as tuples of (type_flag, length).
        self.__cigar = Cigar.__reduce(cigar_sequence)

    @staticmethod
    def from_string(cigar_string):
        raise NotImplementedError()

    def __str__(self):
        return "".join(Cigar.__item_to_string(*item) for item in self.__cigar)

    def __repr__(self):
        return "<Cigar: {!s}>".format(self)

    @staticmethod
    def __item_to_string(flag, length):
        return "{}{}".format(length, Cigar.__get_flag_character(flag))

    @staticmethod
    def __get_flag_character(flag):
        return {
            Cigar.MATCH: "M",
            Cigar.INSERTION: "I",
            Cigar.DELETION: "D",
            Cigar.SKIP: "N",
            Cigar.SOFT_CLIP: "S",
            Cigar.HARD_CLIP: "H",
            Cigar.PADDING: "P",
            Cigar.SEQUENCE_MATCH: "=",
            Cigar.SEQUENCE_MISMATCH: "X"
        }[flag]

    def __add__(self, other):
        assert(isinstance(other, Cigar))
        return Cigar(Cigar.__reduce(self.__cigar + other.__cigar))

    @staticmethod
    def join(*cigars):
        return Cigar(
            Cigar.__reduce(
                itertools.chain(
                    *
                    tuple(
                        item.__cigar for item in cigars))))

    @staticmethod
    def __reduce(cigar_sequence):
        result = []
        iterator = iter(cigar_sequence)
        try:
            prev_item = next(iterator)
        except StopIteration:
            # no items in iterator
            pass
        else:
            for next_item in iterator:
                if prev_item[0] is next_item[0]:
                    # combine
                    prev_item = (prev_item[0], prev_item[1] + next_item[1])
                else:
                    # append
                    result.append(prev_item)
                    prev_item = next_item
            result.append(prev_item)

        # Cigar items like 0M are technically correct but ugly.
        result = [item for item in result if item[1] > 0]
        return result
