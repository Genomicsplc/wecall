# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestAllelePlusStrandBiasFilteringBehaviour(BaseTest):

    def test_should_allow_mildly_strand_biased_calls(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        reads = 10
        strand_bias = 6

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................",
            n_rev=reads + strand_bias, n_fwd=reads - strand_bias, chrom=chrom
        ).with_read(
            "................G...........................",
            n_rev=reads - strand_bias, n_fwd=reads + strand_bias, chrom=chrom
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_no_filters()

    def test_should_allow_mildly_allele_biased_calls(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        reads = 10
        allele_bias = 5

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................",
            n_rev=reads + allele_bias, n_fwd=reads + allele_bias, chrom=chrom
        ).with_read(
            "................G...........................",
            n_rev=reads - allele_bias, n_fwd=reads - allele_bias, chrom=chrom
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_no_filters()

    def test_should_stop_mildly_allele_and_strand_biased_calls(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        reads = 10
        allele_bias = 5
        strand_bias = 4

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            chrom=chrom).with_read(
            "............................................",
            n_rev=reads + allele_bias + strand_bias,
            n_fwd=reads + allele_bias - strand_bias,
            chrom=chrom).with_read(
            "................G...........................",
            n_rev=reads - allele_bias - strand_bias,
            n_fwd=reads - allele_bias + strand_bias,
            chrom=chrom)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_filters({'AB+SB'})

    def test_should_allow_mildly_allele_and_strand_biased_calls_with_lower_specified_threshold(self):
        chrom = 'chr1'
        svc = SVCDriver(self)
        reads = 10
        allele_bias = 5
        strand_bias = 4

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            chrom=chrom).with_read(
            "............................................",
            n_rev=reads + allele_bias + strand_bias,
            n_fwd=reads + allele_bias - strand_bias,
            chrom=chrom).with_read(
            "................G...........................",
            n_rev=reads - allele_bias - strand_bias,
            n_fwd=reads - allele_bias + strand_bias,
            chrom=chrom)
        svc.with_allele_plus_strand_bias_p(0.03)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_no_filters()
