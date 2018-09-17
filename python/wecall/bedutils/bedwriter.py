# All content Copyright (C) 2018 Genomics plc
from wecall.utils.tabix_indexer import TabixIndexer


def bed_line_from_chrom_interval(chrom_interval):
    return "\t".join(
        str(item) for item in [
            chrom_interval.chrom,
            chrom_interval.start,
            chrom_interval.end
        ]
    )


class BEDWriter(object):

    def __init__(self, output_stream):
        self.__output_stream = output_stream

    def write_chrom_interval(self, chrom_interval):
        self.__output_stream.write("{}\n".format(
            bed_line_from_chrom_interval(chrom_interval)))

    def write_bed_record(self, bed_record):
        self.__output_stream.write("{}\n".format(str(bed_record)))

    def write_chrom_intervals(self, chrom_intervals):
        for chrom_interval in chrom_intervals:
            self.write_chrom_interval(chrom_interval)


class BEDWriterContextManager(object):

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        self.fp = open(self.filename, 'w')
        self.vcf_writer = BEDWriter(self.fp)
        return self.vcf_writer

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fp.close()


class BEDIndexer(TabixIndexer):

    def __init__(self, filename):
        TabixIndexer.__init__(self, filename, "bed")
