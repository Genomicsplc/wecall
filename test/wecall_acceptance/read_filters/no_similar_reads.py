# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestNoSimilarReadsFilterFilter(BaseTest):

    def test_should_not_be_on_by_default_in_tests(self):
        svc = SVCDriver(self)
        svc.with_min_reads_per_var(20)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................", n_rev=10, n_fwd=10
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))

    def test_should_filter_similar_reads(self):
        svc = SVCDriver(self)
        svc.with_no_similar_reads_filter(True)
        svc.with_min_reads_per_var(2)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................", n_rev=10, n_fwd=0
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(0)

    def test_should_not_filter_similar_reads_with_different_start_positions(self):
        svc = SVCDriver(self)
        svc.with_no_similar_reads_filter(True)
        svc.with_min_reads_per_var(8)
        svc.with_allow_MNP_calls(False)

        svc.with_ref_sequence(
            # 1234567890123456789
            "TTAATGCATGCATGCATGCATGCATGCATGCATGCCCCG",
        ).with_read(
            "   G...G...G...G...G...G...            ", n_fwd=1, n_rev=1
        ).with_read(
            "       G...G...G...G...G...G...        ", n_fwd=1, n_rev=1
        ).with_read(
            "           G...G...G...G...G...G...    ", n_fwd=1, n_rev=1
        ).with_read(
            ".......G...G...G...G...G...G...........", n_fwd=1, n_rev=1
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(4)

    def test_should_not_filter_similar_reads_with_different_sequences(self):
        svc = SVCDriver(self)
        svc.with_min_reads_per_var(6)
        svc.with_no_similar_reads_filter(True)
        svc.with_allow_MNP_calls(False)

        svc.with_ref_sequence(
            # 1234567890123456789
            "TTAATGCATGCATGCATGCATGCATGCATGCATGCCCCG",
        ).with_read(
            "       ....G...G...G...........        ", n_fwd=1, n_rev=1
        ).with_read(
            "       G.......G...G...........        ", n_fwd=1, n_rev=1
        ).with_read(
            "       G...G.......G...........        ", n_fwd=1, n_rev=1
        ).with_read(
            ".......G...G...G...G...................", n_fwd=1, n_rev=1
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(4)
