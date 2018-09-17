# All content Copyright (C) 2018 Genomics plc
from wecall.genomics import variant
from wecall.vcfutils.record import Record
from wecall.vcfutils.info_data import InfoData
from wecall.vcfutils.sample_data import SampleData
from wecall.vcfutils.schema import Schema
from wecall.vcfutils.writer import VCFWriterContextManager
from wecall_test_drivers.tool_runner import log_file
from wecall.utils.tabix_indexer import TabixIndexer


class VCFBuilder(object):
    def __init__(self, filename, schema=None):
        self.__filename = filename
        self.__indexer = TabixIndexer(self.__filename, "vcf")
        if schema is None:
            self.schema = Schema()
        else:
            self.schema = schema
        self.__records = []

    @property
    def filename(self):
        return self.__filename

    @property
    def compressed_filename(self):
        return self.__indexer.compressed_filename

    @property
    def compressed_filename_index(self):
        return self.__indexer.compressed_filename_index

    def with_variant(self, chrom, pos_from, ref, alt):
        return self.with_record(
            self.generate_record_from_variant(
                variant.Variant(
                    chrom, pos_from, ref, alt)))

    def with_record(self, record):
        self.__records.append(record)
        return self

    def with_record_from_variant(self, variant, **kwargs):
        return self.with_record(
            self.generate_record_from_variant(
                variant, **kwargs))

    def build(self):
        # use default schema
        with VCFWriterContextManager(self.filename, self.schema) as vcf_writer:
            vcf_writer.write_records(self.__records)

        log_file(self.filename)
        return self

    def bgzip(self):
        self.__indexer.bgzip()

    def index(self):
        self.__indexer.index()
        return self

    def generate_record_from_variant(self, variant, **kwargs):
        annotations = {'variant_id': set(),
                       'quality': None,
                       'filters': set(),
                       'info': InfoData(self.schema,
                                        {}),
                       'sample_info': SampleData([key for key,
                                                  _ in self.schema.iter_sample_data()],
                                                 self.schema.samples),
                       'from_multi_alt': False,
                       }
        for key, value in kwargs.items():
            annotations[key] = value

        return Record(schema=self.schema, variant=variant, **annotations)
