# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestAlleleBiasPValue(BaseTest):
    def test_should_get_unknown_value_for_homozygous_alt_call(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT", chrom="1"
        ).with_read(
            ".....................T....................", chrom="1"
        ).with_read(
            ".....................T....................", chrom="1"
        )

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_record_for_variant(Variant("1", 21, "A", "T")) \
            .with_info().with_field("ABPV", [None])

    def test_should_get_unknown_value_for_het_call_with_majority_support(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT", chrom="1"
        ).with_read(
            "..........................................", chrom="1"
        ).with_read(
            ".....................T....................", chrom="1"
        ).with_read(
            ".....................T....................", chrom="1"
        )

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_record_for_variant(Variant("1", 21, "A", "T")) \
            .with_info().with_field("ABPV", [None])


class TestAlleleBiasFilterWithDefaultThreshold(AsciiWecallRunnerTest):

    def test_does_not_filter_homozygous_snp(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............T.............       ",
             "  ...................T........            ",
             "    .................T..............      ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "T", {"DP": [8], "AD": [0, 8]}, ["PASS"])]
        )

    def test_does_not_filter_balanced_heterozygous_snp(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ....................................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "T", {"DP": [8], "AD": [4, 4]}, ["PASS"])]
        )

    def test_does_not_filter_2_out_of_8_heterozygous_snp(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ....................................  ",
             "    ....................................  ",
             "    ....................................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "T", {"DP": [8], "AD": [6, 2]}, ["PASS"])]
        )

    def test_does_not_filter_2_out_of_10_heterozygous_snp(self):
        chrom = '14'
        sample = 'SAMPLE'
        driver = SVCDriver(self)

        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "    ................................      ", n_fwd=4, n_rev=4, chrom=chrom, sample_name=sample
        ).with_read(
            "    .................T..................  ", n_fwd=1, n_rev=1, chrom=chrom, sample_name=sample)

        vcf = driver.call().with_output_vcf()

        vcf \
            .has_record_for_variant(Variant(chrom, 21, 'A', 'T')) \
            .with_filters(set()) \
            .with_sample(sample) \
            .has_read_depth(10) \
            .has_allelic_read_support(8, 2)

    def test_does_not_filter_2_out_of_12_heterozygous_snp(self):
        chrom = '14'
        sample = 'SAMPLE'
        driver = SVCDriver(self)

        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "    ................................      ", n_fwd=5, n_rev=5, chrom=chrom, sample_name=sample
        ).with_read(
            "    .................T..................  ", n_fwd=1, n_rev=1, chrom=chrom, sample_name=sample)

        vcf = driver.call().with_output_vcf()

        vcf \
            .has_record_for_variant(Variant(chrom, 21, 'A', 'T')) \
            .with_filters(set()) \
            .with_sample(sample) \
            .has_read_depth(12) \
            .has_allelic_read_support(10, 2)

    def test_does_filter_2_out_of_15_heterozygous_snp(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAGTACACATACGCACGCGCCAGCACGTGAATTGATCTTGTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ....................................  ",
             "    ....................................  ",
             "    ....................................  ",
             "    .................T..................  ",
             ".....................T......              "],
            []  # We don't call this because the variant coverage is too low
        )

    def test_does_not_filter_5_out_of_40_heterozygous_snp(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "T", {"DP": [40], "AD": [35, 5]}, ["AB"])],
            config_dict={"varFilterIDs": "AB"}
        )


class TestAlleleBiasFilterWithNonDefaultThresholds(AsciiWecallRunnerTest):

    def test_filters_2_out_ot_12_heterozygous_snp(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ....................................  ",
             "    ....................................  ",
             "    ....................................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "T", {"DP": [12], "AD": [10, 2]}, ["AB"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_does_not_filter_3_out_ot_12_heterozygous_snp(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ....................................  ",
             "    ....................................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "T", {"DP": [12], "AD": [9, 3]}, ["PASS"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_filters_2_out_ot_12_heterozygous_insertion(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCC**TAAAAAAAATTTTTTTTTTT",
            ["       .............**.............       ",
             "  ..................**........            ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..................  ",
             "    ................**..................  ",
             "    ................**..................  ",
             "    ................TT..................  ",
             "....................TT......              "],
            [(19, "C", "CTT", {"DP": [12], "AD": [10, 2]}, ["AB"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_does_not_filter_3_out_ot_12_heterozygous_insertion(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCC**TAAAAAAAATTTTTTTTTTT",
            ["       .............**.............       ",
             "  ..................**........            ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..............      ",
             "    ................**..................  ",
             "    ................**..................  ",
             "    ................TT..................  ",
             "    ................TT..................  ",
             "....................TT......              "],
            [(19, "C", "CTT", {"DP": [12], "AD": [9, 3]}, ["PASS"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_filters_2_out_ot_12_heterozygous_deletion(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCGTTAAAAAAAATTTTTTTTTTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ....................................  ",
             "    ....................................  ",
             "    ....................................  ",
             "    ................**..................  ",
             "....................**......              "],
            [(19, "CGT", "C", {"DP": [12], "AD": [10, 2]}, ["AB"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_does_not_filter_3_out_ot_12_heterozygous_deletion(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCCGTTAAAAAAAATTTTTTTTTTT",
            ["       ............................       ",
             "  ............................            ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ................................      ",
             "    ....................................  ",
             "    ....................................  ",
             "    ................**..................  ",
             "    ................**..................  ",
             "....................**......              "],
            [(19, "CGT", "C", {"DP": [12], "AD": [9, 3]}, ["PASS"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_filters_3_out_ot_18_heterozygous_insertion_in_A4_homopolymer(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCC**AAAAGCCGTTTTTTTTTTT",
            ["       .............**............       ",
             "  ..................**...............    ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.................  ",
             "    ................**.................  ",
             "    ................**.................  ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.................  ",
             "    ................**.................  ",
             "    ................**.................  ",
             "    ................AA.................  ",
             "    ................AA.................  ",
             "....................AA.....              "],
            [(19, "C", "CAA", {"DP": [18], "AD": [15, 3]}, ["AB"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_does_not_filter_3_out_ot_12_heterozygous_insertion_in_A4_homopolymer(self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCC**AAAAGCCGTTTTTTTTTTT",
            ["       .............**............       ",
             "  ..................**.......            ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.............      ",
             "    ................**.................  ",
             "    ................**.................  ",
             "    ................AA.................  ",
             "    ................AA.................  ",
             "....................AA.....              "],
            [(19, "C", "CAA", {"DP": [12], "AD": [9, 3]}, ["PASS"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )

    def test_does_not_filter_3_out_ot_12_heterozygous_insertion_in_A4_homopolymer_when_some_variants_need_left_aligning(
            self):
        self.calls_variants_with_sample_data_and_filters(
            "AAAAAAAAAAACGCACCCCC**AAAA**GCCGTTTTTTTTTTT",
            ["       .............**....**........       ",
             "  ..................**....**...            ",
             "    ................**....**.........      ",
             "    ................**....**.........      ",
             "    ................**....**.........      ",
             "    ................**....**.........      ",
             "    ................**....**.........      ",
             "    ................**....**.............  ",
             "    ................**....**.............  ",
             "    ................**....AA.............  ",
             "    ................AA....**.............  ",
             "....................AA....                 "],
            [(19, "C", "CAA", {"DP": [12], "AD": [9, 3]}, ["PASS"])],
            config_dict={"minAlleleBiasP": "0.04"}
        )
