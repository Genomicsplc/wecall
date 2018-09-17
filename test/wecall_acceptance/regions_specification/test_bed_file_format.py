# All content Copyright (C) 2018 Genomics plc
import os
from unittest import expectedFailure

from wecall.bedutils.bedwriter import BEDWriterContextManager
from wecall.utils.interval import ChromInterval
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestBedFileInput(BaseTest):

    @expectedFailure
    def test_bed_extension_is_not_required(self):
        region_filename = os.path.join(self.work_dir, "bed.txt")

        with BEDWriterContextManager(region_filename) as bed:
            bed.write_chrom_interval(ChromInterval("1", 0, 29))

        svc = SVCDriver(self).with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT", chrom="1"
        ).with_read(
            "....G...................", n_rev=10, n_fwd=10
        ).with_region_string(region_filename)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1)

    def test_should_run_if_interval_not_contained_in_reference(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT", chrom="1"
        ).with_read(
            "....G...................", chrom="1"
        ).with_bed_file(
            ['1\t24\t28']
        ).with_region_padding(0)

        svc.call(True).with_output_vcf().record_count(0)

    def test_should_warn_user_if_contigs_provided_not_in_reference(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT", chrom="1"
        ).with_read(
            "....G...................", chrom="1"
        ).with_bed_file(["42\t1\t5"])

        expect = svc.call(True)

        expect.with_log().bed_file_contains_contigs_that_are_not_present_in_the_reference_warning("42")

    def test_should_warn_even_if_some_contigs_not_in_reference(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT", chrom="1"
        ).with_read(
            "....G...................", chrom="1"
        ).with_bed_file(["1\t0\t10", "42\t21881\t22032"])

        expect = svc.call(True)

        expect.with_log().bed_file_contains_contigs_that_are_not_present_in_the_reference_warning("42")

    def test_ok_if_multiple_regions_not_contained_in_reference(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT", chrom="1"
        ).with_read(
            "....G...................", chrom="1"
        ).with_bed_file(
            ["21\t0\t10", "42\t21881\t22032", "42\t22032\t22033"]
        )

        expect = svc.call(True)

        expect.with_log().bed_file_contains_contigs_that_are_not_present_in_the_reference_warning("21", "42")

    def test_disallow_mixing_bed_files_and_region_strings(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT"
        ).with_read(
            "....G..................."
        ).with_region_string("some_bed_file.bed,1:1-10")

        expect = svc.call(False)

        expect.regions_contains_both_bedfile_and_region_string_error()

    def test_disallow_mixing_bed_files_and_region_strings2(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT"
        ).with_read(
            "....G..................."
        ).with_region_string("1:1-10,some_bed_file.bed")

        expect = svc.call(False)

        expect.regions_contains_both_bedfile_and_region_string_error()

    def test_error_if_input_regions_file_does_not_exist(self):
        bed_file_name = "some_bed_file.bed"

        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "ACGTACGTACGTACGTACGTACGT"
        ).with_read(
            "....G..................."
        ).with_region_string(bed_file_name)

        expect = svc.call(False)

        expect.bedfile_does_not_exist_error(bed_file_name)
