# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestLowQualityFilter(BaseTest):
    def test_should_filter_low_quality_call(self):
        chrom = 'chr1'
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "................G...........................",
            n_rev=1, n_fwd=1, chrom=chrom
        ).with_read(
            "............................................",
            n_rev=1, n_fwd=1, chrom=chrom
        ).with_min_call_qual(40)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_filters({'LQ'})

    def test_should_not_filter_high_quality_call(self):
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
        ).with_min_call_qual(40)

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(chrom, 16, 'T', 'G')) \
            .with_no_filters()
