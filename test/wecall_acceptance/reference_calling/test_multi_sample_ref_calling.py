# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestReferenceCallingMultiSample(BaseTest):
    def test_calls_whole_region_as_reference_if_no_variants(self):
        chrom = "1"
        sample_1 = "sample_1"
        sample_2 = "sample_2"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ...........................       ", chrom=chrom, sample_name=sample_1, n_fwd=5, n_rev=5
        ).with_read(
            "  ................................       ", chrom=chrom, sample_name=sample_2, n_fwd=5, n_rev=5
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 41)

    def test_calls_correct_reference_when_one_sample_has_snp_and_tother_has_indel(self):
        chrom = "1"
        sample_1 = "sample_1"
        sample_2 = "sample_2"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACG*CCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..........T................        ", chrom=chrom, sample_name=sample_1
        ).with_read(
            "  ...............T...........             ", chrom=chrom, sample_name=sample_1
        ).with_read(
            "    ...........T.*................        ", chrom=chrom, sample_name=sample_2
        ).with_read(
            "    ...........T.*.....................   ", chrom=chrom, sample_name=sample_2
        ).with_read(
            "...............T.*.........               ", chrom=chrom, sample_name=sample_2
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()

        # Has only 4 records which are:-
        vcf_expect.has_reference_calls_for_region(chrom, 0, 15)
        vcf_expect.has_record(chrom, 15, "C", "T")
        vcf_expect.has_record(chrom, 16, "G", "GT")
        vcf_expect.has_reference_calls_for_region(chrom, 17, 41)

    def test_calls_correct_ref_calls_with_one_snp_even_if_reference_for_tother(self):
        chrom = "1"
        sample_1 = "sample_1"
        sample_2 = "sample_2"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..........T................       ", chrom=chrom, sample_name=sample_1
        ).with_read(
            "  ...............T...........            ", chrom=chrom, sample_name=sample_1
        ).with_read(
            "    ...............................      ", chrom=chrom, sample_name=sample_2
        ).with_read(
            "    ...................................  ", chrom=chrom, sample_name=sample_2
        ).with_read(
            "...........................              ", chrom=chrom, sample_name=sample_2
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 17)
        vcf_expect.has_record(chrom, 17, "C", "T").with_sample(sample_2).has_genotype("0/0")
        vcf_expect.has_reference_calls_for_region(chrom, 18, 41)
