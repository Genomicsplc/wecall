# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.read_sequence import FORWARD_GOOD_READ
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestPairedReadFilters(BaseTest):

    def test_should_call_variant_with_proper_paired_reads_and_allow_improper_reads_flag_set(self):
        svc = SVCDriver(self)
        svc.with_allow_improper_pairs()

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

    def test_should_call_variant_with_improper_paired_reads_and_allow_improper_reads_flag_set(self):
        svc = SVCDriver(self)
        svc.with_allow_improper_pairs()

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................", n_rev=10, n_fwd=10, read_flags=FORWARD_GOOD_READ & ~2
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1) \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))

    def test_should_call_variant_with_improper_paired_reads_when_allow_improper_reads_flag_not_set(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
        ).with_read(
            "................G...........................", n_rev=10, n_fwd=10, read_flags=FORWARD_GOOD_READ & ~2
        )

        expect = svc.call()

        expect.with_output_vcf() \
            .record_count(1)
