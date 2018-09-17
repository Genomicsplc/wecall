# All content Copyright (C) 2018 Genomics plc
import re
import functools

from wecall.genomics.chromosome import chromosome_comp
from wecall.utils.interval import ChromInterval

"""
The variant module defines the standardised, internal representation of variant data which
is used throughout the echidna code.
"""


# Module constants
TYPE_REF = 0
TYPE_SNP = 1
TYPE_INS = 2
TYPE_DEL = 3
TYPE_MNP = 4
TYPE_SYM = 5

SMALL_VARIANT_TYPES = [TYPE_REF, TYPE_SNP, TYPE_INS, TYPE_DEL, TYPE_MNP]

TYPE_TO_STR = {
    TYPE_REF: 'REF',
    TYPE_SNP: 'SNP',
    TYPE_INS: 'INS',
    TYPE_DEL: 'DEL',
    TYPE_MNP: 'MNP',
    TYPE_SYM: 'SYM'
}

STR_TO_TYPE = {value: key for key, value in list(TYPE_TO_STR.items())}

DB_VARIANT_TYPES = {
    "All": {
        TYPE_SNP,
        TYPE_INS,
        TYPE_DEL,
        TYPE_MNP
    },
    "SNP": {TYPE_SNP},
    "MNP": {TYPE_MNP},
    "Ins": {TYPE_INS},
    "Del": {TYPE_DEL},
    "INDEL": {TYPE_DEL, TYPE_INS},
}


def mnp_to_snps(var):
    """
    Split multi-nucleotide replacement variant objects into multiple
    objects representing SNPs.

    Apparently this is deprectated - for what reason and what should it be replaced with?
    """
    assert var.type == TYPE_MNP, "Invalid variant type"

    for index, (refBase, altBase) in enumerate(zip(var.ref, var.alt)):
        if refBase != altBase:
            yield Variant(var.chrom, var.pos_from + index, refBase, altBase)


def variant_type(ref, alt):
    """
    Determine the type of variant we have, based on the REF and ALT
    strings:
    """
    # all REFs
    if ref == alt or alt == '.' or alt == '<NON_REF>':
        return TYPE_REF

    # Symbolic alleles - (typically ignored or filtered out)
    if alt.startswith(('<')):
        return TYPE_SYM

    # simple SNPs
    if len(ref) == 1 and len(alt) == 1:
        return TYPE_SNP

    # INDELs
    elif len(ref) < len(alt):
        return TYPE_INS
    elif len(ref) > len(alt):
        return TYPE_DEL

    # MNPs or complex SNPs
    else:
        # NOTE: len(ref) == len(alt) here.
        n_diff = len([(c1, c2) for c1, c2 in zip(ref, alt) if c1 != c2])
        assert(n_diff > 0)
        return TYPE_SNP if n_diff == 1 else TYPE_MNP


@functools.total_ordering
class Variant(object):
    """
    Standard internal variant representation. To be used everywhere in the
    echidna code when we are storing and passing around variant data.
    """

    __slots__ = ('chrom', 'pos_from', 'ref', 'alt')

    def __init__(self, chrom, pos_from, ref, alt):
        self.chrom = chrom
        self.pos_from = pos_from
        self.ref = ref
        self.alt = alt

    @property
    def pos_to(self):
        return self.pos_from + len(self.ref)

    @property
    def insert_size(self):
        return len(self.alt) - self.length

    @property
    def length(self):
        return len(self.ref)

    @property
    def type(self):
        return variant_type(self.ref, self.alt)

    @property
    def one_indexed_pos_from(self):
        return self.pos_from + 1

    @property
    def one_indexed_pos_to(self):
        return self.pos_to + 1

    def overlap(self, chrom_interval):
        return chrom_interval.overlap(
            ChromInterval(
                self.chrom,
                self.pos_from,
                self.pos_to))

    def __hash__(self):
        """
        This allows the Variant to be used as e.g. a key for a
        dictionary, or stored in a set.
        """
        return hash((
            self.chrom,
            self.pos_from,
            self.pos_to,
            self.ref,
            self.alt
        ))

    def __eq__(self, other):
        """
        NOTE: This is non-symmetric (by def self cannot be None)
        """
        # TODO: deal with ambiguous representations caused by non-left-aligned
        # data
        if other is None:
            return False

        return self.chrom == other.chrom and\
            self.pos_from == other.pos_from and\
            self.pos_to == other.pos_to and\
            self.ref == other.ref and\
            self.alt == other.alt

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        """
        Define comparison operator for variants, to allow sorting.
        """
        # TODO: deal with ambiguous representations caused by non-left-aligned
        # data
        if self.chrom != other.chrom:
            return chromosome_comp(self.chrom, other.chrom)
        elif self.pos_from != other.pos_from:
            return self.pos_from < other.pos_from
        elif self.pos_to != other.pos_to:
            return self.pos_to < other.pos_to
        elif self.ref != other.ref:
            return self.ref < other.ref
        elif self.alt != other.alt:
            return self.alt < other.alt
        else:
            return False

    def __str__(self):
        return "Variant({}: {} - {}, {} --> {})".format(self.chrom,
                                                        self.pos_from, self.pos_to, self.ref, self.alt)

    def __repr__(self):
        return self.__str__()
