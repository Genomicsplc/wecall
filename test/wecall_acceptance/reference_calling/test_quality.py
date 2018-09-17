# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver
from wecall_test_drivers.vcf_expectation import ref_alt


class TestRefCallingQuality(BaseTest):
    def test_get_unknown_quality_if_no_reads_span_region(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "                                         ",
        ).with_output_ref_calls(True)

        driver.call().with_output_vcf().has_record(chrom, 0, "A", ref_alt).with_quality(None)
