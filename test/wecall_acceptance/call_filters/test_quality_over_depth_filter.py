# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver

MAX_PHRED = 3000


class TestQualityOverDepthFilter(BaseTest):

    def test_should_apply_soft_filter_to_snp_with_low_quality(self):
        svc = SVCDriver(self)
        svc.with_var_filters("QD")
        svc.with_min_snp_q_over_depth(MAX_PHRED)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................", n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        record_expectation = vcf_expectation.has_record_for_variant(
            Variant(DEFAULT_CHROM, 16, "T", "G"))
        record_expectation.with_filters({"QD"})

    def test_should_not_apply_soft_filter_to_snp_with_high_quality(self):
        svc = SVCDriver(self)
        svc.with_var_filters("QD")
        svc.with_min_snp_q_over_depth(1.0)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................", n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        record_expectation = vcf_expectation.has_record_for_variant(
            Variant(DEFAULT_CHROM, 16, "T", "G"))
        record_expectation.with_filters(set())

    def test_should_apply_soft_filter_to_indel_with_low_quality(self):
        svc = SVCDriver(self)
        svc.with_var_filters("QD")
        svc.with_min_indel_q_over_depth(MAX_PHRED)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGG*TAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................C...........................",
            n_rev=1, n_fwd=1
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        record_expectation = vcf_expectation.has_record_for_variant(
            Variant(DEFAULT_CHROM, 15, "G", "GC"))
        record_expectation.with_filters({"QD"})

    def test_should_not_apply_soft_filter_to_indel_with_high_quality(self):
        svc = SVCDriver(self)
        svc.with_var_filters("QD")
        svc.with_min_indel_q_over_depth(1.0)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGG*TAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................C...........................", n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        record_expectation = vcf_expectation.has_record_for_variant(
            Variant(DEFAULT_CHROM, 15, "G", "GC"))
        record_expectation.with_filters(set())

    def test_should_not_apply_filter_for_quality_below_cap(self):
        record_expectation = QD_impl(self, n_rev=30, n_fwd=30)
        record_expectation.with_filters(set())

    def test_should_not_apply_filter_for_quality_beyond_cap(self):
        # Quality is capped so that quality / (number of supporting reads) is
        # low artificially.
        record_expectation = QD_impl(self, n_rev=50, n_fwd=50)
        record_expectation \
            .with_quality(MAX_PHRED) \
            .with_filters(set())
        record_expectation.with_info() \
            .with_field('QD', [None])


def QD_impl(test_case, n_fwd, n_rev):
    svc = SVCDriver(test_case)
    svc.with_var_filters("QD")
    svc.with_min_snp_q_over_depth(35)

    svc.with_ref_sequence(
        # 1234567890123456789
        "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
    ).with_read(
        "................G...........................", n_rev=n_rev, n_fwd=n_fwd
    )

    expect = svc.call()
    vcf_expectation = expect.with_output_vcf()
    vcf_expectation.record_count(1)

    return vcf_expectation.has_record_for_variant(
        Variant(DEFAULT_CHROM, 16, "T", "G"))
