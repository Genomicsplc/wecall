# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.read_sequence import FORWARD_GOOD_READ, DUPLICATE
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestNoDuplicatesFilter(BaseTest):

    def test_should_call_variant_if_reads_are_not_duplicates(self):
        svc = SVCDriver(self)
        svc.with_duplicates_filter(True)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................",
            n_rev=10, n_fwd=10, read_flags=FORWARD_GOOD_READ & ~DUPLICATE
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))

    def test_should_not_call_variant_with_duplicate_reads_and_duplicate_reads_not_allowed(self):
        svc = SVCDriver(self)
        svc.with_duplicates_filter(True)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................",
            n_rev=10, n_fwd=10, read_flags=FORWARD_GOOD_READ | DUPLICATE
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(0)

    def test_should_not_call_variant_with_duplicate_reads_with_default(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................",
            n_rev=10, n_fwd=10, read_flags=FORWARD_GOOD_READ | DUPLICATE
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(0)

    def test_should_call_variant_with_duplicate_reads_and_duplicate_reads_allowed(self):
        svc = SVCDriver(self)
        svc.with_duplicates_filter(False)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................",
            n_rev=10, n_fwd=10, read_flags=FORWARD_GOOD_READ | DUPLICATE
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))
