# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


MAX_PHRED = 3000


class TestMultiSampleDiploidPhasedCalls(BaseTest):
    def test_phasing_for_isolated_snp_on_one_sample_only(self):
        sample_1 = "sample_1"
        sample_2 = "sample_2"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
            "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom
        ).with_read(
            "........T............", n_fwd=10, n_rev=10, sample_name=sample_1
        ).with_read(
            ".....................", n_fwd=10, n_rev=10, sample_name=sample_2
        )
        svc_driver.with_output_phased_genotypes(True)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(1)

        record_expect = vcf_expect.has_record_for_variant(
            Variant(chrom, 8, "G", "T"))

        sample_1_expect = record_expect.with_sample(sample_1)
        sample_1_expect.has_phased_genotype("1|1")
        sample_1_expect.has_phase_set_id(str(8))
        sample_1_expect.has_phase_set_quality(MAX_PHRED)

        sample_2_expect = record_expect.with_sample(sample_2)
        sample_2_expect.has_phased_genotype("0|0")
        sample_2_expect.has_phase_set_id(str(8))
        sample_2_expect.has_phase_set_quality(MAX_PHRED)

    def test_phasing_for_two_snps_across_samples(self):
        sample_1 = "sample_1"
        sample_2 = "sample_2"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
            "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom
        ).with_read(
            "........T............", n_fwd=10, n_rev=10, sample_name=sample_1
        ).with_read(
            ".....................", n_fwd=10, n_rev=10, sample_name=sample_1
        ).with_read(
            "........C............", n_fwd=10, n_rev=10, sample_name=sample_2
        ).with_read(
            ".....................", n_fwd=10, n_rev=10, sample_name=sample_2
        )
        svc_driver.with_output_phased_genotypes(True)

        expect = svc_driver.call()\
            .with_output_vcf() \
            .record_count(2)

        expect.has_record_for_variants(Variant(chrom, 8, "G", "T"), Variant(chrom, 8, "G", "C"))
        expected_phase_set_id = str(8)

        expect \
            .has_record_for_variants(Variant(chrom, 8, "G", "T"), Variant(chrom, 8, "G", "C")) \
            .with_sample(sample_1)\
            .has_phased_genotypes("0|1", "0|.") \
            .has_phase_set_id(expected_phase_set_id) \
            .has_phase_set_quality(MAX_PHRED)

        expect \
            .has_record_for_variants(Variant(chrom, 8, "G", "T"), Variant(chrom, 8, "G", "C")) \
            .with_sample(sample_2)\
            .has_phased_genotypes(".|0", "1|0") \
            .has_phase_set_id(expected_phase_set_id) \
            .has_phase_set_quality(MAX_PHRED)

    def test_phasing_for_overlapping_variant_across_samples(self):
        sample_1 = "sample_1"
        sample_2 = "sample_2"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
            "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom
        ).with_read(
            "........T............", n_fwd=10, n_rev=10, sample_name=sample_1
        ).with_read(
            "........T............", n_fwd=10, n_rev=10, sample_name=sample_2
        ).with_read(
            "........C............", n_fwd=10, n_rev=10, sample_name=sample_2
        )
        svc_driver.with_output_phased_genotypes(True)

        expect = svc_driver.call()\
            .with_output_vcf() \
            .record_count(2)

        expected_phase_set_id = str(8)

        expect \
            .has_record_for_variants(Variant(chrom, 8, "G", "C"), Variant(chrom, 8, "G", "T")) \
            .with_sample(sample_1)\
            .has_phased_genotypes(".|.", "1|1") \
            .has_phase_set_id(expected_phase_set_id) \
            .has_phase_set_quality(MAX_PHRED)

        expect \
            .has_record_for_variants(Variant(chrom, 8, "G", "C"), Variant(chrom, 8, "G", "T")) \
            .with_sample(sample_2)\
            .has_phased_genotypes("1|.", ".|1") \
            .has_phase_set_id(expected_phase_set_id) \
            .has_phase_set_quality(MAX_PHRED)
