# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestDefaultStrandBiasFilteringBehaviour(BaseTest):

    def test_should_allow_unbiased_calls_to_pass_through(self):
        chrom = 'chr1'
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "................G...........................",
            n_rev=10, n_fwd=10, chrom=chrom
        ).with_read(
            "............................................",
            n_rev=10, n_fwd=10, chrom=chrom
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_no_filters()

    def test_should_stop_forward_biased_calls_to_pass_through(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        reads = 10
        bias = 7

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................",
            n_rev=reads + bias, n_fwd=reads - bias, chrom=chrom
        ).with_read(
            "................G...........................",
            n_rev=reads - bias, n_fwd=reads + bias, chrom=chrom
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_filters({'SB'})

    def test_should_stop_reverse_biased_calls_to_pass_through(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        reads = 10
        bias = -7

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................",
            n_rev=reads + bias, n_fwd=reads - bias, chrom=chrom
        ).with_read(
            "................G...........................",
            n_rev=reads - bias, n_fwd=reads + bias, chrom=chrom
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_filters({'SB'})


class TestCustomStrandBiasThreshold(BaseTest):
    """
    Estimating custom probability cut-off point
    > def beta_binomial(k, n, a, b):
    .     kCn = scipy.special.comb(n, k)
    .     num = scipy.special.beta(k+a, n-k+b)
    .     den = scipy.special.beta(a, b)
    .     return kCn * num / den
    > sum((beta_binomial(n, 30, 20.0*30.0/30.0, 20.0) for n in range(21, 30)))
    0.062450575207941214
    > sum((beta_binomial(n, 30, 20.0*30.0/30.0, 20.0) for n in range(22, 30)))
    0.033780564203978708
    > 0.5 * (0.062450575207941214 + 0.033780564203978708)
    0.04811556970595996
    """

    def test_should_allow_forward_biased_calls_just_below_custom_threshold_to_pass_through(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        svc.with_strand_bias_p(0.04811556970595996)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................", n_rev=21, n_fwd=9, chrom=chrom
        ).with_read(
            "................G...........................", n_rev=9, n_fwd=21, chrom=chrom)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_no_filters()

    def test_should_filter_forward_biased_calls_just_above_custom_threshold(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        svc.with_strand_bias_p(0.04811556970595996)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................", n_rev=22, n_fwd=8, chrom=chrom
        ).with_read(
            "................G...........................", n_rev=8, n_fwd=22, chrom=chrom)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_filters({'SB'})

    def test_should_allow_reverse_biased_calls_just_below_custom_threshold_to_pass_through(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        svc.with_strand_bias_p(0.04811556970595996)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................", n_rev=9, n_fwd=21, chrom=chrom
        ).with_read(
            "................G...........................", n_rev=21, n_fwd=9, chrom=chrom)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_no_filters()

    def test_should_filter_reverse_biased_calls_just_above_custom_threshold(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        svc.with_strand_bias_p(0.04811556970595996)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................", n_rev=22, n_fwd=8, chrom=chrom
        ).with_read(
            "................G...........................", n_rev=8, n_fwd=22, chrom=chrom)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_filters({'SB'})
