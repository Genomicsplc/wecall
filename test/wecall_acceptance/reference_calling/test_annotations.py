# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


ref_alt = "<NON_REF>"


class TestRefCallingMinDepthComputation(BaseTest):
    def test_depth_computation_all_reads_spanning_reference(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            ".........................................", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_output_ref_calls(True)

        expect = driver.call()

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

    def test_depth_computation_all_reads_spanning_reference_with_one_snp(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "................T........................", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_output_ref_calls(True).with_allow_MNP_calls(False)

        expect = driver.call()
        vcf_expect = expect.with_output_vcf()

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 17, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

    def test_depth_computation_all_reads_spanning_reference_with_variants(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "................T.T......................", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_output_ref_calls(True).with_allow_MNP_calls(False)

        expect = driver.call()
        vcf_expect = expect.with_output_vcf()

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 17, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 19, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

    def test_depth_computation_all_reads_spanning_reference_with_insertion(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAC*AAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "................T.......................", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_output_ref_calls(True).with_allow_MNP_calls(False)

        expect = driver.call()
        vcf_expect = expect.with_output_vcf()

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 15, "C", "CT"))\
            .with_sample(sample_name).has_read_depth(10)

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 16, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

    def test_depth_computation_all_reads_spanning_reference_with_deletion(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAACACAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "................*.......................", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_output_ref_calls(True).with_allow_MNP_calls(False)

        expect = driver.call()
        vcf_expect = expect.with_output_vcf()

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 15, "CA", "C"))\
            .with_sample(sample_name).has_read_depth(10)

        vcf_expect \
            .has_record_for_variant(Variant(chrom, 17, "C", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

    def test_min_depth_computation_with_mixed_depth_of_reads(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "..............................          ", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_read(
            "          ..............................", n_rev=3, n_fwd=3, sample_name=sample_name
        ).with_output_ref_calls(True)

        expect = driver.call()

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 10, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(16).has_min_read_depth(16)

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 30, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(6).has_min_read_depth(6)

    def test_min_depth_computation_start_boundary_conditions(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            " .......................................", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_output_ref_calls(True)

        expect = driver.call()

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(0).has_min_read_depth(0)

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 1, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

    def test_min_depth_computation_end_boundary_conditions(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "....................................... ", n_rev=5, n_fwd=5, sample_name=sample_name
        ).with_output_ref_calls(True)

        expect = driver.call()

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(10).has_min_read_depth(10)

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 39, "A", ref_alt))\
            .with_sample(sample_name).has_read_depth(0).has_min_read_depth(0)

    def test_min_depth_computation_with_mixed_depth_of_reads_when_no_chunking_occurs(self):
        sample_name = "bah"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", chrom=chrom
        ).with_read(
            "........................................", n_rev=10, n_fwd=10, sample_name=sample_name
        ).with_read(
            "...................................     ", n_rev=2, n_fwd=2, sample_name=sample_name
        ).with_read(
            "   ...................................  ", n_rev=1, n_fwd=1, sample_name=sample_name
        ).with_output_ref_calls(True)

        expect = driver.call()

        expect\
            .with_output_vcf()\
            .has_record_for_variant(Variant(chrom, 0, "A", ref_alt))\
            .with_sample(sample_name)\
            .has_read_depth(round(20 + 4 * 35 / 40 + 2 * 35 / 40))\
            .has_min_read_depth(20)
