# All content Copyright (C) 2018 Genomics plc
from wecall.bedutils.bedwriter import BEDWriterContextManager
from wecall.utils.interval import ChromInterval
from wecall.wecall_utils.log_utils import log_timing_parser
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver
import os


class TestVariantCallerTimings(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.bed_filename = os.path.join(self.work_dir, "_.bed")
        self.log_filename = os.path.join(self.work_dir, "_.log")

        chrom = "1"

        with BEDWriterContextManager(self.bed_filename) as bed_file:
            bed_file.write_chrom_interval(ChromInterval(chrom, 0, 42))

        self.output_vcf = os.path.join(self.work_dir, "output.vcf")
        self.svc_driver = SVCDriver(self)
        self.svc_driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT", chrom=chrom
        )
        self.svc_driver.with_read(
            "       ..............C.............       ", chrom=chrom
        )
        self.svc_driver.with_log_timings(True)
        self.svc_driver.with_region_string(self.bed_filename)
        self.svc_driver.with_log_filename(self.log_filename)
        self.svc_driver.with_output_vcf_filename(self.output_vcf)

    def __outputs_timings_for_files(self, expected_files):
        observed_files = set()

        with open(self.log_filename, "r") as log_file:
            timing_data = log_timing_parser(log_file)
            self.assertGreater(len(timing_data), 0)
            for timing_data_item in timing_data:
                self.assertEqual("IO", timing_data_item.timing_type)
                self.assertEqual("us", timing_data_item.length_units)
                self.assertIn("file", timing_data_item.metadata)
                observed_files.add(timing_data_item.metadata["file"])

        for expected in expected_files:
            self.assertIn(expected, observed_files)

    def test_should_contain_timings_output_for_bam(self):
        self.svc_driver.with_bam_filenames(
            [os.path.join(self.work_dir, "ba.bam")])
        self.svc_driver.call()
        self.__outputs_timings_for_files(
            {os.path.join(self.work_dir, "ba.bam")})

    def test_should_contain_timings_output_for_fasta_file(self):
        ref = os.path.join(self.work_dir, "ref.fa")
        self.svc_driver.with_ref_filename(ref)

        self.svc_driver.call()
        self.__outputs_timings_for_files({ref})

    def test_should_contain_timings_output_for_fasta_index_file(self):
        ref = os.path.join(self.work_dir, "ref.fa")
        self.svc_driver.with_ref_filename(ref)

        self.svc_driver.call()
        self.__outputs_timings_for_files({ref + ".fai"})

    def test_should_contain_timings_output_for_vcf(self):
        self.svc_driver.call()
        self.__outputs_timings_for_files({self.output_vcf})

    def test_should_contain_timings_output_for_vcf_when_variant_caller_run_in_parallel(self):
        self.svc_driver.with_number_of_jobs(1)
        self.svc_driver.with_work_dir(
            os.path.join(self.work_dir, "vc_work_dir"))

        self.svc_driver.call()
        self.__outputs_timings_for_files({self.output_vcf})

    def test_should_contain_timings_output_for_bed(self):
        self.svc_driver.call()
        self.__outputs_timings_for_files({self.bed_filename})
