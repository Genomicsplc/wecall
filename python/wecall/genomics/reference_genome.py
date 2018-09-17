# All content Copyright (C) 2018 Genomics plc
from abc import ABCMeta, abstractmethod
from wecall.genomics.chromosome import CHROMOSOME_ORDER
from wecall.genomics.reference_chromosome import ReferenceChromosome


class AbstractReferenceGenome(object, metaclass=ABCMeta):

    @abstractmethod
    def chromosomes(self):
        return []

    @abstractmethod
    def fetch(self, chrom, start, end):
        return ''

    @abstractmethod
    def get_chrom_length(self, chrom):
        return 0


class InMemoryReferenceGenome(AbstractReferenceGenome):

    def __init__(self):
        self.__data = {}

    def with_chrom(self, name, sequence, pos_from=0):
        reference_chrom = self.__data[name] = ReferenceChromosome(
            sequence, pos_from, name)
        return reference_chrom

    def chromosomes(self):
        return sorted(list(self.__data.keys()),
                      key=lambda x: CHROMOSOME_ORDER.get(x) or -1)

    def get_chrom_length(self, chrom):
        return self.__data[chrom].pos_to

    def fetch(self, chrom, start=None, end=None):
        if start is None:
            start = 0
        if end is None:
            end = self.get_chrom_length(chrom)

        seq = self.__data[chrom].fasta_string()[start:end]
        if end - start != len(seq):
            raise IndexError
        return seq
