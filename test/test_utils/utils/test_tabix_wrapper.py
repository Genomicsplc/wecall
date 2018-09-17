# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
import wecall.utils.interval
from wecall.vcfutils.writer import VCFWriterContextManager
from os import path, remove
from wecall.utils.tabix_wrapper import TabixWrapper
from wecall_test_drivers.base_test import BaseTest
from wecall.utils.tabix_indexer import TabixIndexer


class TestTabixWrapper(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.vcf_output = path.join(self.work_dir, self.id() + ".vcf")
        with VCFWriterContextManager(self.vcf_output) as vcf_writer:
            vcf_writer.write_variant(Variant("20", 61097, "C", "T"))

        indexer = TabixIndexer(self.vcf_output, "vcf")
        indexer.index()
        self.vcf_output = indexer.compressed_filename
        self.tabix_file = TabixWrapper(self.vcf_output)

    def tearDown(self):
        remove(self.vcf_output)
        remove(self.vcf_output + ".tbi")

    def test_should_yield_records_from_tabix_file_with_standard_chrom_interval(self):
        records = list(
            self.tabix_file.fetch_generator(
                wecall.utils.interval.ChromInterval(
                    "20", 61097, 61098)))
        self.assertEqual(len(records), 1)

    def test_should_not_yield_any_records_for_region_after_record(self):
        records = list(
            self.tabix_file.fetch_generator(
                wecall.utils.interval.ChromInterval(
                    "20", 61098, 70000)))
        self.assertEqual(len(records), 0)

    def test_should_not_yield_any_records_for_region_before_record(self):
        records = list(
            self.tabix_file.fetch_generator(
                wecall.utils.interval.ChromInterval(
                    "20", 60000, 61097)))
        self.assertEqual(len(records), 0)

    def test_should_yield_regions_for_whole_chromosome_chrom_interval_format(self):
        records = list(self.tabix_file.fetch_generator(
            wecall.utils.interval.ChromInterval("20")))
        self.assertEqual(len(records), 1)

    def test_should_raise_stop_iteration_for_invalid_chromosome(self):
        chrom_interval = wecall.utils.interval.ChromInterval("Hello")

        try:
            next(self.tabix_file.fetch_generator(chrom_interval))
            raised = False
        except StopIteration:
            raised = True

        self.assertTrue(raised)
