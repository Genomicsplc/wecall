# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver
from math import log10


MAX_PHRED = 3000


class TestSingleSampleDiploidPhasedCalls(BaseTest):
    def test_phasing_for_isolated_homozygous_alt_variant(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom
        ).with_read(
            "........T............", n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_output_phased_genotypes(True)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(1)

        record_expect = vcf_expect.has_record_for_variant(Variant(chrom, 8, "G", "T"))

        sample_expect = record_expect.with_sample(sample_name)
        sample_expect.has_genotype("1|1")
        sample_expect.has_phase_set_id(str(8))
        sample_expect.has_phase_set_quality(MAX_PHRED)

    def test_phasing_for_isolated_heterozygous_variant(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom) \
            .with_read(
                "........T............", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".....................", n_fwd=10, n_rev=10, sample_name=sample_name)

        svc_driver.with_output_phased_genotypes(True)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(1)

        record_expect = vcf_expect.has_record_for_variants(
            Variant(chrom, 8, "G", "T"))

        sample_expect = record_expect.with_sample(sample_name)
        sample_expect.has_phased_genotypes("0|1")
        sample_expect.has_phase_set_id(str(8))
        sample_expect.has_phase_set_quality(MAX_PHRED)

    def test_phasing_for_two_heterozygous_variants_ocrn_same_strand(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self) \
            .with_ref_sequence(
                "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom) \
            .with_read(
                "........T...T........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".....................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        records_expect = vcf_expect.has_record_for_variants(
            Variant(chrom, 8, "G", "T"),
            Variant(chrom, 12, "A", "T")
        )
        records_expect\
            .with_sample(sample_name)\
            .has_phased_genotypes("0|1", "0|1")\
            .has_phase_set_id("8")

    def test_phasing_for_two_heterozygous_variants_on_same_strand_more_than_a_cluster_apart(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom) \
            .with_read(
                "........TT...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".....................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True)\
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(-1)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect .has_record_for_variants(
            Variant(
                chrom, 8, "G", "T"), Variant(
                chrom, 9, "C", "T")) .with_sample(sample_name) .has_phased_genotypes(
            "0|1", "0|1")

    def test_phasing_for_two_heterozygous_variants_on_different_strand(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom) \
            .with_read(
                "........T............", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "............T........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True)\
            .with_allow_MNP_calls(False)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        records_expect = vcf_expect.has_record_for_variants(
            Variant(chrom, 8, "G", "T"),
            Variant(chrom, 12, "A", "T")
        )
        records_expect.with_sample(sample_name)\
            .has_phased_genotypes("0|1", "1|0")\
            .has_phase_set_id("8")

    @expectedFailure
    def test_symmetry_in_repetitive_reference(self):
        sample_name = "a_sample"
        chrom = "1"

        m = 2
        svc_driver = SVCDriver(self) \
            .with_ref_sequence(
                "TAAAAAAAAAAAAAAAAAAAAAAAAAT", chrom=chrom) \
            .with_read(
                "........T..................", n_fwd=m, n_rev=m, sample_name=sample_name) \
            .with_read(
                "...........................", n_fwd=m, n_rev=m, sample_name=sample_name) \
            .with_read(
                ".................T.........", n_fwd=m, n_rev=m, sample_name=sample_name) \
            .with_allow_MNP_calls(False)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect.has_record_for_variants(Variant(chrom, 17, "A", "T"))\
            .with_sample(sample_name)\
            .has_phased_genotypes("0/1")

        vcf_expect.has_record_for_variants(Variant(chrom, 8, "A", "T"))\
            .with_sample(sample_name)\
            .has_phased_genotypes("0/1")

    def test_phase_quality_for_equally_likely_phasings_for_same_variant_representations(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom) \
            .with_read(
                "........T...T..T.....", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "............T........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "........T...T........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "............T..T.....", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True)\
            .with_allow_MNP_calls(False)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(3)

        # for het calls ratio_of_phase_to_total is defined as
        # N_supporting_mutation / ( N_supporting_mutation + N_supporting_ref)
        ratio_of_phase_to_total = 0.5
        expected_phase_quality = int(
            round(log10(1.0 - ratio_of_phase_to_total) * -10.0))

        vcf_expect.has_record_for_variants(Variant(chrom, 8, "G", "T"))\
            .with_sample(sample_name)\
            .has_phased_genotypes("0|1")\
            .has_phase_set_id("8")\
            .has_phase_set_quality(expected_phase_quality)

        vcf_expect.has_record_for_variants(Variant(chrom, 12, "A", "T"))\
            .with_sample(sample_name)\
            .has_phased_genotypes("1|1")\
            .has_phase_set_id("8")\
            .has_phase_set_quality(expected_phase_quality)

        vcf_expect.has_record_for_variants(Variant(chrom, 15, "A", "T"))\
            .with_sample(sample_name)\
            .has_phased_genotypes("0|1")\
            .has_phase_set_id("8")\
            .has_phase_set_quality(expected_phase_quality)

    @expectedFailure
    def test_phase_quality_for_phase_with_2_out_of_3_support(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCTGCAAAAAAAAAAT", chrom=chrom)\
            .with_read(
                "........T...T..T.....", n_fwd=5, n_rev=5, sample_name=sample_name) \
            .with_read(
                "............T........", n_fwd=5, n_rev=5, sample_name=sample_name) \
            .with_read(
                "........T...T........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "............T..T.....", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True)\
            .with_allow_MNP_calls(False)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(3)

        ratio_of_phase_to_total = 2.0 / 3.0
        # actual value needs to be figured out from equations
        unknown_phase_quality = int(round(log10(1.0 - ratio_of_phase_to_total) * -10.0))

        vcf_expect.has_record_for_variants(
            Variant(chrom, 8, "G", "T"),
            Variant(chrom, 12, "A", "T"),
            Variant(chrom, 15, "A", "T")
        )\
            .with_sample(sample_name)\
            .has_phased_genotypes("0|1", "1|1", "1|0")\
            .has_phase_set_id("8")\
            .has_phase_set_quality(unknown_phase_quality)
