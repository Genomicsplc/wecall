# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall_test_drivers.ascii_quality_recalibration_runner import AsciiQualityRecalibrationTest
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest


class TestQualityRecalibrationDeletion(AsciiQualityRecalibrationTest):
    def test_should_not_recalibrate_good_read_data_for_deletion(self):
        reference = "ATCTAATAGCTATCAGCAATATCGCGCGTATTATTTATTTAT"
        bam_spec = ["    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,     ",
                    "..............*.............              ",
                    "      ,,,,,,,,*,,,,,,,,,,,,,,,,,,,,       ",
                    "..............*.......................    "]

        self.assert_quality_recalibrated_in_output_bam(reference, bam_spec, bam_spec)


class TestQualityRecalibrationInsertion(AsciiQualityRecalibrationTest):
    def test_should_not_recalibrate_good_read_data_for_deletion(self):
        reference = "ATCTAATAGCTATC*GCAATATCGCGCGTATTATTTATTTAT"
        bam_spec = ["    ,,,,,,,,,,a,,,,,,,,,,,,,,,,,,,,,,     ",
                    "..............A.............              ",
                    "      ,,,,,,,,a,,,,,,,,,,,,,,,,,,,,       ",
                    "..............A.......................    "]

        self.assert_quality_recalibrated_in_output_bam(reference, bam_spec, bam_spec)


class TestQualityRecalibrationBespokeQualities(AsciiQualityRecalibrationTest):
    def test_should_not_recalibrate_good_read_data_with_snp_1(self):
        reference = "ATCTAATAGCTATCAGCAATATCGCGCGTATTATTTATTTAT"
        bam_spec = ["    ,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,     ",
                    "                            0             ",
                    "..............T.............              ",
                    "                           0              ",
                    "      ,,,,,,,,t,,,,,,,,,,,,,,,,,,,,       ",
                    "       0                                  ",
                    "..............T.......................    ",
                    " 0                                        "]

        self.assert_quality_recalibrated_in_output_bam(reference, bam_spec, bam_spec)

    def test_should_not_recalibrate_good_read_data_with_snp_2(self):
        reference = "ATCTAATAGCTATCAGCAATATCGCGCGTATTATTTATTTAT"
        bam_spec = ["    ,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,     ",
                    "                            1             ",
                    "..............T.............              ",
                    "                           1              ",
                    "      ,,,,,,,,t,,,,,,,,,,,,,,,,,,,,       ",
                    "       1                                  ",
                    "..............T.......................    ",
                    " 1                                        "]

        self.assert_quality_recalibrated_in_output_bam(reference, bam_spec, bam_spec)


class TestQualityRecalibrationSingleSNP(AsciiQualityRecalibrationTest):
    def test_should_not_recalibrate_region_snp_on_two_forward_strands_out_of_eight(self):
        reference = "ATCTAATAGCATCTAATAGCTAGCATCCGTAACAGCAATATCGCGCGTATTATTTATTTAT"
        bam_spec = [
            "..............................T..............................",
            "         .....................T..............................",
            "         ....................................................",
            "         ....................................................",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            "     ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "]

        self.assert_quality_recalibrated_in_output_bam(reference, bam_spec, bam_spec)

    def test_should_not_recalibrate_snp_on_one_forward_and_reverse_strand_out_of_eight(self):
        reference = "ATCTAATAGCATCTAATAGCTAGCATCCGTAACAGCAATATCGCGCGTATTATTTATTTAT"
        bam_spec = [
            "..............................T..........................    ",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            "     ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "]

        self.assert_quality_recalibrated_in_output_bam(reference, bam_spec, bam_spec)

    def test_should_recalibrate_around_snp_on_two_forward_strands_out_of_nine(self):
        reference = "ATCTAATAGCATCTAATAGCTAGCATCCGTAACAGCAATATCGCGCGTATTATTTATTTAT"
        input_bam = [
            "..............................T..............................",
            "         .....................T..............................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            "     ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "]

        output_bam = [
            "..............................T..............................",
            "                   000000000000000000000000000000000000000000",
            "         .....................T..............................",
            "                     0000000000000000000000000000000000000000",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            "     ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "]

        self.assert_quality_recalibrated_in_output_bam(reference, input_bam, output_bam)

    def test_should_recalibrate_snp_on_one_forward_and_reverse_strand_out_of_nine(self):
        reference = "ATCTAATAGCATCTAATAGCTAGCATCCGTAACAGCAATATCGCGCGTATTATTTATTTAT"
        input_bam = [
            "..............................T..............................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            "     ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "]

        output_bam = [
            "..............................T..............................",
            "                   000000000000000000000000000000000000000000",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            "         ....................................................",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ",
            " 0000000000000000000000000000000000000000                    ",
            "     ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "]

        self.assert_quality_recalibrated_in_output_bam(reference, input_bam, output_bam)


class TestRealDataExamplesFromNA12878(AsciiWecallRunnerTest):

    def calls_variants_without_recalibration(self, ref, sequence_list, expected_ascii_haplotypes):
        self.calls_variants(
            ref, sequence_list,
            config_dict={"recalibrateBaseQs": False, "overwrite": True},
            expected_ascii_haplotypes=expected_ascii_haplotypes
        )

    def calls_variants_with_recalibration(self, ref, sequence_list, expected_ascii_haplotypes):
        self.calls_variants(
            ref, sequence_list,
            config_dict={"recalibrateBaseQs": True, "overwrite": True},
            expected_ascii_haplotypes=expected_ascii_haplotypes
        )

    def calls_variants_with_and_without_recalibration(self, ref, sequence_list, expected_ascii_haplotypes):
        self.calls_variants_with_recalibration(
            ref, sequence_list, expected_ascii_haplotypes)
        self.calls_variants_without_recalibration(
            ref, sequence_list, expected_ascii_haplotypes)

    def test_calls_two_good_snps_with_and_without_recalibration(self):
        self.calls_variants_with_and_without_recalibration(
            "ACGCCCCCTGCAAAAACTACTAAAAA",
            [".T........................",
             ".T........................",
             "...........C..............",
             "...........C.............."],
            [".T........................",  # Expected calls
             "...........C.............."]
        )

    @expectedFailure
    def test_calls_false_positive_snp_with_and_without_recalibration(self):
        self.calls_variants_without_recalibration(
             "AGTGCCTGTTGCAAACTTAAAGTAT**********AA**********TAAAATAAA**********ATAAATAAAAAAAAATAAAAAAAAGAATA",
            [",,,,,,,,,,,    ..........**********..**********.........**********.............................",
             ".................  ......**********..**********.........**********.............................",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,   ...**********.............................",
             ".........................**********..**********.....  ..**********.............................",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,                 ,,,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,      ..................",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,,,,    .................",
             "...............T.........**********G.**********.........**********..............    ...........",
             "               1                   1                                                           ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,  ..........",
             ",,,t,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,        ................",
             "   1                                                                                           ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,aataaaataa,,**********,,,,,,,,,**********,,,,,,,,,,,,,        ........",
             "                         3333333333                                                            ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,                                         ,,,,,",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,,,,,,,,,,,,,,,,,,aataaaataa,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,  ",
             "                         3333333333                                                            ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,aataaaataa,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,  ",
             "                         3333333333                                                            ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             ",,aa,,,,,,,,,,,,,,,,,a,,,aataaaataa,,**********,,,,,,,,,**********,,,,,g,,,,,,,,,,,,,,,,,,,,,, ",
             "  11                 1   3333333333                                    1                       ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,taaaataaac,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "                                     3333333333                                                ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,aataaaataa,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "                         3333333333                                                            ",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,**********,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,,,,,,,,,,,,,,,,,,**********,,**********,,,,,,,,,ataaaataaa,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "                                                        3333333333                             ",
             ".........................**********..**********.........**********.....AT....T............   ,,",
             "                                                                             1                 ",
             ".........................**********..**********.........**********.............................",
             ".........................**********..**********.........**********.....AT....T................."],

            [".........................**********..**********.........**********.....AT......................",  # Expected calls  # noqa
             ".........................AATAAAATAA..**********.........**********.....AT......................"]  # Expected calls  # noqa
        )
