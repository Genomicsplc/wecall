# All content Copyright (C) 2018 Genomics plc
from wecall.utils.interval import ChromInterval
import re
from wecall.common.exceptions import EchidnaException
from wecall.genomics.chromosome import standardise_chromosome


DEFAULT_CHROM = "1"


class ReferenceChromosome(object):

    def __init__(self, ref_string, pos_from=0, chrom=DEFAULT_CHROM):
        self.__validate_ref_seq(ref_string)

        self.chrom = chrom
        self.pos_from = pos_from
        self.ref_seq = ref_string
        self._ref_minus_deletions = self.ref_seq.replace("*", "")

    def __str__(self):
        return self._ref_minus_deletions

    def length_with_deletions(self):
        return len(self.ref_seq)

    def length_minus_deletions(self):
        return len(self._ref_minus_deletions)

    def __getitem__(self, item):
        return self._ref_minus_deletions[item - self.pos_from]

    @property
    def chrom_interval(self):
        return ChromInterval(self.chrom, self.pos_from, self.pos_to)

    @property
    def pos_to(self):
        return self.pos_from + len(self._ref_minus_deletions)

    def fasta_string(self):
        return self.pos_from * 'N' + self._ref_minus_deletions

    def __validate_ref_seq(self, ref_seq):
        if not re.match(r'^[ACGTURYKMSWBDHVN\*]*\Z', ref_seq):
            raise EchidnaException(
                "Illegal character in reference sequence {!r}".format(ref_seq))
