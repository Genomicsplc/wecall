# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestVariantClustering(AsciiWecallRunnerTest):

    def calls_variants_with_coverage_20(
            self, ref, sequence_list, expected_variants=None):
        self.calls_variants(
            ref, sequence_list, expected_ascii_haplotypes=None,
            expected_variant_stubs=expected_variants,
            n_fwd=20, n_rev=20
        )

    def test_calls_eight_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAAATT",
            [".T.*.A.*.A..*.G.*................."]
        )

    def test_calls_nine_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAATT",
            [".T.*.A.*.A..*.G.*.C.............."]
        )

    def test_calls_ten_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAATT",
            [".T.*.A.*.A..*.G.*.C.*..........."]
        )

    def test_calls_eleven_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T........"]
        )

    def test_calls_twelve_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*......"]
        )

    def test_calls_thirteen_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAAAAAAAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T........."]
        )

    def test_calls_fourteen_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAATAAAAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T.*......."]
        )

    def test_calls_fifteen_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAATAAAAAAAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T.*.T........"]
        )

    def test_calls_sixteen_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAATAAATAACAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T.*.T.*......"]
        )

    def test_calls_seventeen_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAATAAATAACAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T.*.T.*.T...."]
        )

    def test_calls_eighteen_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAATAAATAAATAACAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T.*.T.*.T.*......."]
        )

    def test_calls_nineteen_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAATAAATAAATAAACAAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T.*.T.*.T.*.T......."]
        )

    def test_calls_twenty_variants(self):
        self.calls_variants_with_coverage_20(
            "ACGCTCACTGCAGTCGTTAGAAATAAATAAATAAATAAATACAAATT",
            [".T.*.A.*.A..*.G.*.C.*.T.*.T.*.T.*.T.*.T.*......"]
        )


class TestSNPClustering(AsciiWecallRunnerTest):

    def calls_variants_with_coverage_20_no_MNPs(
            self, ref, sequence_list, expected_variants=None):
        self.calls_variants(
            ref, sequence_list,
            config_dict={"allowMNPCalls": False},
            expected_ascii_haplotypes=None,
            expected_variant_stubs=expected_variants,
            n_fwd=20, n_rev=20
        )

    def test_calls_eight_SNPs(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAAAAAAAA",
            [".T.T.A.G.A..C.G.A...................."]
        )

    def test_calls_eight_adjecent_SNPs(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAAAAAAAA",
            [".TCTAAGGC............................"]
        )

    def test_calls_nine_SNPs(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAAAAAAAA",
            [".T.T.A.G.A..C.G.A..C................."]
        )

    def test_calls_nine_adjacent_SNPs(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAAAAAAAA",
            ["............CGGTACGCT................"]
        )

    def test_calls_twenty_adjacent_SNPs(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGTTAGTGCAATGAAAAAAAAAA",
            ["........ACTCCGTACGCTACGTGCTC........."]
        )


class TestCallingTwoNearbyClusters(AsciiWecallRunnerTest):

    def calls_variants_with_coverage_20_no_MNPs(
            self, ref, sequence_list, expected_variants=None):
        self.calls_variants(
            ref, sequence_list,
            config_dict={"maxClusterDist": 5, "allowMNPCalls": False},
            expected_ascii_haplotypes=None,
            expected_variant_stubs=expected_variants,
            n_fwd=20, n_rev=20
        )

    def test_calls_two_sets_of_eight_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTTCCGTCGTTAGAAAAAAAAAAAAAAAAA",
            [".TGTGATGG.............GCTTCGTC............"]
        )

    def test_calls_two_sets_of_nine_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTTCCGTCGTTAGAAAAAAAAAAAAAAAAA",
            [".TGTGATGGT............GCTTCGTCT..........."]
        )

    def test_calls_two_sets_of_ten_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTTCCGTCGTTAGAAAAAAAAAAAAAAAAA",
            [".TGTGATGGTG...........GCTTCGTCTG.........."]
        )

    def test_calls_two_sets_of_eleven_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCCGTAGTTAGAAAAAAAAAAAAAAAAA",
            [".TGTGATGGTGT.........GCTTCGTCTGT........."]
        )

    def test_calls_two_sets_of_twelve_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGCTGGTTAGAAAAAAAAAAAAAAAAA",
            [".TGTGATGGTGTC........GCTTCGTCTGTC........"]
        )


class TestCallingThreeNearbyClusters(AsciiWecallRunnerTest):

    def calls_variants_with_coverage_20_no_MNPs(
            self, ref, sequence_list, expected_variants=None):
        self.calls_variants(
            ref, sequence_list,
            config_dict={"maxClusterDist": 5, "allowMNPCalls": False},
            expected_ascii_haplotypes=None,
            expected_variant_stubs=expected_variants,
            n_fwd=20, n_rev=20
        )

    def test_calls_three_sets_of_eight_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGCTGACTGATCGCTGGTTAGAAAAAAAACCTGACTGAGTACAAAAAAAAACGTAGTACGTA",
            [".TGTGATGG....................GCTTCGTC..................GGTCTCGT.........."])

    def test_calls_three_sets_of_nine_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGCCTGACTGATGGTTAGAAAAAAAACGCTGACTGATACAAAAAAAAACGTAGTACGTA",
            [".TGTGATGGT...................GCTTCGTCT.................GGTCTCGTT........."])

    def test_calls_three_sets_of_ten_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "AAAAAAAACGCTCACTGCAGTCGCCTGACTGATGGTTAGAAAAAAAACGTCTGACTGAACAAAAAAAAACGTAGTACGTA",
            ["       .TGTGATGGTG..................GCTTCGTCTG................GGTCTCGTTG........"])

    def test_calls_three_sets_of_eleven_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGCCTGACTGATGGTTAGAAAAAAAACGTACCTGACTGAAAAAAAAAACGTAGTACGTA",
            [".TGTGATGGTGG.................GCTTCGTCTGG...............GGTCTCGTTGG......."])

    def test_calls_three_sets_of_twelve_SNPs_in_separate_clusters(self):
        self.calls_variants_with_coverage_20_no_MNPs(
            "ACGCTCACTGCAGTCGCTCTGACTGAGGTTAGAAAAAAAACGTACTGACTGACAAAAAAAAACGTAGTACGTA",
            [".TGTGATGGTGGT................GCTTCGTCTGGT..............GGTCTCGTTGGT......"])


class TestCallingOnVeryLongVariantClusters(BaseTest):

    def test_call_on_ridiculous_snp_cluster(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AGCTAGCGCTAGCGCTCGACAGATCGAGATAGCCGGGCTAAGATTAGATCGCGATGCGATGCACGTACGCATGCATACGA"
        ).with_read(
            "CAGGCATATGCATATGTACTCACGTACACGCATTAAATGCCACGGCACGTATACGATACGATCTAGCTATCGATCGCTAC", n_fwd=20, n_rev=20
        ).with_allow_MNP_calls(False)

        expect = driver.call()

        expect \
            .with_output_vcf() \
            .record_count(80)


class TestCallingWithIndelsInRepetitiveSequence(BaseTest):
    def test_indel_that_gets_left_aligned_over_snp_gets_correct_genotype(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "ATCGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAATCG"
        ).with_read(
            "............................................T*******....", n_fwd=10, n_rev=10, sample_name="sample"
        ).with_read(
            "..G******...............................................", n_fwd=10, n_rev=10, sample_name="sample")

        expect = driver.call()
        expect\
            .with_output_vcf()\
            .has_record_for_variants(
                Variant("1", 1, "TCGATTA", "T"),
                Variant("1", 2, "CGATTACA", "C"),
                Variant("1", 8, "C", "G"),
                Variant("1", 51, "A", "T"),
            ).with_sample("sample").has_phased_genotypes('1|.', '.|1', '1|.', '0|1')


class TestAllowMNPCallsParameter(BaseTest):

    def test_calls_SNPs_when_not_allowed_to_call_MNPs(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAAAAAAAA"
        ).with_read(
            "............CGGTACGCT................", n_fwd=20, n_rev=20
        ).with_allow_MNP_calls(False)

        expect = driver.call()

        expect \
            .with_output_vcf() \
            .record_count(9)

    def test_calls_MNPs_when_allowed_to_call_MNPs(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "ACGCTCACTGCAGTCGTTAGAAAAAAAAAAAAAAAAA"
        ).with_read(
            "............CGGTACGCT................", n_fwd=20, n_rev=20
        ).with_allow_MNP_calls(True)

        expect = driver.call()

        expect \
            .with_output_vcf() \
            .record_count(1)
