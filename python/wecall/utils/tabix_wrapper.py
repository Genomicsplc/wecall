# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.chromosome import standardise_chromosome
import pysam


class TabixWrapper(object):

    def __init__(self, tabix_filename):
        self.__tabix_file = pysam.Tabixfile(tabix_filename, 'r')
        self.__contig_mapping = {standardise_chromosome(
            contig): contig for contig in self.__tabix_file.contigs}

    @property
    def header(self):
        return (line for line in self.__tabix_file.header)

    @property
    def contigs(self):
        return self.__tabix_file.contigs

    def fetch_generator(self, chrom_interval):
        # Tabix will throw a ValueError if the chromosome specified is not
        # present in the index for this file.
        try:
            if chrom_interval.chrom is None:
                return self.__tabix_file.fetch()
            else:
                return self.__tabix_file.fetch(
                    self.__contig_mapping.get(
                        chrom_interval.chrom,
                        chrom_interval.chrom),
                    chrom_interval.interval.start,
                    chrom_interval.interval.end)
        except ValueError:
            raise StopIteration

    def fetch_region(self, region):
        try:
            return self.__tabix_file.fetch(region=region)
        except ValueError:
            raise StopIteration

    def close(self):
        self.__tabix_file.close()

    def __enter__(self):
        return self

    def __exit__(self, ex_type, value, traceback):
        self.close()
