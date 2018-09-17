# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestMinSquareRootMeanSquareMappingQualityFilter(BaseTest):
    def test_should_filter_variant_when_all_reads_have_quality_below_threshold(self):
        svc = SVCDriver(self)
        threshold = 50
        svc.with_var_filters("MQ")
        svc.with_min_root_mean_square_mapping_q(threshold)
        svc.with_read_mapping_filter_q(0)
        svc\
            .with_ref_sequence(
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG")\
            .with_read(
                "................G...........................", n_fwd=10, n_rev=10, mapping_quality=threshold - 1)

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation\
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))\
            .with_filters({"MQ"})

    def test_should_filter_variant_when_all_supporting_reads_have_low_mapping_quality(self):
        svc = SVCDriver(self)
        threshold = 50
        low_mq = threshold - 1
        high_mq = threshold * 2
        svc.with_var_filters("MQ")
        svc.with_min_root_mean_square_mapping_q(threshold)
        svc.with_read_mapping_filter_q(0)
        svc\
            .with_ref_sequence(
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG")\
            .with_read(
                "................G.........................  ", n_fwd=10, n_rev=10, mapping_quality=low_mq)\
            .with_read(
                "............................................", n_fwd=10, n_rev=10, mapping_quality=high_mq)

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation\
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))\
            .with_filters({"MQ"})

    def test_should_not_filter_variant_when_all_supporting_reads_have_high_mapping_quality(self):
        svc = SVCDriver(self)
        threshold = 50
        low_mq = threshold - 1
        high_mq = threshold
        svc.with_var_filters("MQ")
        svc.with_min_root_mean_square_mapping_q(threshold)
        svc.with_read_mapping_filter_q(0)
        svc\
            .with_ref_sequence("AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG").with_read(
                "................G.........................  ", n_fwd=10, n_rev=10, mapping_quality=high_mq)\
            .with_read(
                "............................................", n_fwd=10, n_rev=10, mapping_quality=low_mq)

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)
        vcf_expectation\
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))\
            .with_no_filters()

    def test_should_take_max_mq_over_all_the_samples_which_support_the_variant(self):
        svc = SVCDriver(self)
        threshold = 50
        svc.with_var_filters("MQ")
        svc.with_min_root_mean_square_mapping_q(threshold)
        svc.with_read_mapping_filter_q(0)
        svc\
            .with_ref_sequence(
                # 1234567890123456789
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG")\
            .with_read(
                "................G.........................  ",
                n_fwd=10, n_rev=10, mapping_quality=threshold - 1, sample_name="Ugly")\
            .with_read(
                "................G...........................",
                n_fwd=10, n_rev=10, mapping_quality=threshold, sample_name="Good")

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)
        vcf_expectation\
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))\
            .with_no_filters()

    def test_should_take_max_mq_only_over_the_the_samples_that_support_the_variant(self):
        svc = SVCDriver(self)
        threshold = 50
        svc.with_var_filters("MQ")
        svc.with_min_root_mean_square_mapping_q(threshold)
        svc.with_read_mapping_filter_q(0)
        svc\
            .with_ref_sequence(
                # 1234567890123456789
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG")\
            .with_read(
                "................G.........................  ",
                n_fwd=10, n_rev=10, mapping_quality=threshold - 1, sample_name="Ugly")\
            .with_read(
                "............................................",
                n_fwd=10, n_rev=10, mapping_quality=threshold + 1, sample_name="Good")

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)
        vcf_expectation\
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))\
            .with_filters({"MQ"})

    def test_should_filter_variant_when_root_mean_square_of_supporting_reads_is_below_threshold(self):
        svc = SVCDriver(self)
        threshold = 50
        svc.with_var_filters("MQ")
        svc.with_min_root_mean_square_mapping_q(threshold)
        svc.with_read_mapping_filter_q(0)
        svc\
            .with_ref_sequence(
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG")\
            .with_read(
                "................G...........................", n_fwd=10, n_rev=10, mapping_quality=threshold - 1)\
            .with_read(
                "................G...........................", n_fwd=10, n_rev=10, mapping_quality=threshold)

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation\
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))\
            .with_filters({"MQ"})

    def test_should_not_filter_variant_when_root_mean_square_of_supporting_reads_is_above_threshold(self):
        svc = SVCDriver(self)
        threshold = 100
        svc.with_var_filters("MQ")
        svc.with_min_root_mean_square_mapping_q(threshold)
        svc.with_read_mapping_filter_q(0)
        svc\
            .with_ref_sequence(
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG")\
            .with_read(
                "................G...........................", n_fwd=10, n_rev=10, mapping_quality=threshold - 1)\
            .with_read(
                "................G...........................", n_fwd=1, n_rev=1, mapping_quality=threshold + 30)

        expect = svc.call()
        vcf_expectation = expect.with_output_vcf()
        vcf_expectation.record_count(1)

        vcf_expectation\
            .has_record_for_variant(Variant(DEFAULT_CHROM, 16, "T", "G"))\
            .with_no_filters()
