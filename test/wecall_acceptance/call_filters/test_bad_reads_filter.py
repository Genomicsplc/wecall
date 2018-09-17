# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestBadReadsFilter(BaseTest):
    def test_bad_reads_filter_applied_when_snp_has_low_quality_bases_on_left_of_supporting_reads(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(7) \
            .with_min_bad_reads_score(21)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
            "         2                                  ",
            n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G")) \
            .with_filters({"BR"})

    def test_bad_reads_filter_applied_when_snp_has_low_quality_bases_on_right_of_supporting_reads(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(7) \
            .with_min_bad_reads_score(21)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
            "                       2                    ",
            n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G")) \
            .with_filters({"BR"})

    def test_bad_reads_filter_window_considers_full_alignment_span_of_indel(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(7) \
            .with_min_bad_reads_score(21)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................*...........................",
            "                        2                   ",
            n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 15, "GT", "G")) \
            .with_filters({"BR"})

    def test_bad_reads_filter_window_considers_full_alignment_span_of_indel_in_di_nucleotide_region(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(7) \
            .with_min_bad_reads_score(21)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTGTGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................**..........................",
            "                           2                ",
            n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 14, "GGT", "G")) \
            .with_filters({"BR"})

    def test_should_not_apply_filter_to_snp_if_all_supporting_reads_are_good(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(7) \
            .with_min_bad_reads_score(15)

        svc.with_ref_sequence(
            # 1234567  890123456789
            "AAAGCGTAA**CCGGGTTAGT**CAAACCCGTTACGTATGCATG"
        ).with_read(
            ".........**.....G....**.....................", n_rev=10, n_fwd=10
        ).with_read(
            ".........GT..........TA.....................",
            "       00               00                  ",
            n_rev=11, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(3)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 14, "T", "G")) \
            .with_no_filters()

    def test_should_not_apply_bad_reads_to_insertion_if_all_supporting_reads_have_high_base_qualities(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(3) \
            .with_min_bad_reads_score(15)

        svc.with_ref_sequence(
            # 1234567890123 456789
            "AAAGCGTACAACCG*GGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "..............*..G...........................",
            "                 1                           ",
            n_rev=11, n_fwd=10
        )
        svc.with_read(
            "..............T..............................",
            n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 13, "G", "GT")) \
            .with_no_filters()

    def test_should_not_apply_filter_with_base_pair_too_far_on_left_of_snp(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(7) \
            .with_min_bad_reads_score(15)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
            "        0                                   ", n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G")) \
            .with_no_filters()

    def test_should_not_apply_filter_with_base_pair_too_far_on_right_of_snp(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(1) \
            .with_min_bad_reads_score(15)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
            "                  0                         ", n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G")) \
            .with_no_filters()

    def test_bad_reads_filter_not_applied_when_median_read_is_good(self):
        svc = SVCDriver(self) \
            .with_var_filters("BR") \
            .with_bad_reads_window_size(7) \
            .with_min_bad_reads_score(20)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
            "           1      1                         ", n_rev=10, n_fwd=10
        ).with_read(
            "................G...........................",
            "         4444444 4444444                    ", n_rev=11, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G")) \
            .with_no_filters()

    def test_bad_reads_filter_not_applied_when_snp_has_high_quality_bases_nearby(self):
        svc = SVCDriver(self)
        svc.with_var_filters("BR")
        svc.with_bad_reads_window_size(7)
        svc.with_min_bad_reads_score(15)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
            "         4444444 4444444                    ",
            n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G")) \
            .with_no_filters()

    def test_bad_reads_filter_not_applied_if_one_sample_is_not_naughty(self):
        svc = SVCDriver(self)
        svc.with_var_filters("BR")
        svc.with_bad_reads_window_size(7)
        svc.with_min_bad_reads_score(13)

        svc.with_ref_sequence(
            # 1234567890123456789
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG"
        ).with_read(
            "................G...........................",
            "         3333333 3333333                    ",
            sample_name="GOOD", n_rev=2, n_fwd=2
        ).with_read(
            "................G...........................",
            "         0000000 0000000                    ",
            sample_name="BAD", n_rev=10, n_fwd=10
        ).with_read(
            "................G...........................",
            "         00000      0000                    ",
            sample_name="UGLY", n_rev=10, n_fwd=10
        )

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation \
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G")) \
            .with_no_filters()
