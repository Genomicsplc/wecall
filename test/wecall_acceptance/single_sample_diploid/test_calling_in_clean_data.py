# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestCallingWithDifferentReadLengths(AsciiWecallRunnerTest):

    def test_should_not_call_snps_with_reads_of_length_less_or_equal_to_20(self):
        self.calls_variants(
            "ACGCCCCTGCAAAAAAAAAA",
            [".T..................",
             ".T..................",
             "..........C.........",
             "..........C........."],
            ["....................",
             "...................."]
        )

    def test_calls_two_snps_on_reads_longer_than_20(self):
        self.calls_variants(
            "ACGCCCCCTGCAAAAAAAAAA",
            [".T...................",
             ".T...................",
             "...........C.........",
             "...........C........."]
        )


class TestCallingSingleVariants(AsciiWecallRunnerTest):

    def test_calls_heterozygous_snp_on_the_left_edge(self):
        self.calls_variants_with_genotype(
            "ACGCCCCCTGCAAAAAAAAAA",
            ["C....................",
             "C....................",
             ".....................",
             "....................."],

            ["C....................",
             "....................."]  # Expected genotype
        )

    def test_calls_homozygous_snp_on_the_left_edge(self):
        self.calls_variants_with_genotype(
            "ACGCCCCCTGCAAAAAAAAAA",
            ["C....................",
             "C....................",
             "C....................",
             "C...................."],

            ["C....................",
             "C...................."]  # Expected genotype
        )

    def test_calls_heterozygous_snp_on_the_right_edge(self):
        self.calls_variants_with_genotype(
            "ACGCCCCCTGCAATAAATGCG",
            ["....................C",
             ".....................",
             ".....................",
             "....................C"],

            ["....................C",
             "....................."]  # Expected genotype
        )

    def test_calls_homozygous_snp_on_the_right_edge(self):
        self.calls_variants_with_genotype(
            "ACGCCCCCTGCGACAGATAAT",
            ["....................C",
             "....................C",
             "....................C",
             "....................C"],

            ["....................C",
             "....................C"]  # Expected genotype
        )

    def test_calls_heterozygous_snp(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            ["...................T........................",
             "...................T........................",
             "............................................",
             "............................................"],

            ["...................T........................",  # Expected genotype
             "............................................"]
        )

    def test_calls_homozygous_snp(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            ["...................T........................",
             "...................T........................",
             "...................T........................",
             "...................T........................"],

            ["...................T........................",  # Expected genotype
             "...................T........................"]
        )

    def test_calls_heterozygous_insertion(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTA*GTCACAAACCCGTTACGTATGCATG",
            ["...................T.........................",
             "...................T.........................",
             "...................*.........................",
             "...................*........................."],

            ["...................T.........................",  # Expected genotype
             "...................*........................."]
        )

    def test_calls_homozygous_insertion(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTA*GTCACAAACCCGTTACGTATGCATG",
            ["...................T.........................",
             "...................T.........................",
             "...................T.........................",
             "...................T........................."],

            ["...................T.........................",  # Expected genotype
             "...................T........................."]
        )

    def test_calls_heterozygous_deletion(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            ["............................................",
             "............................................",
             "...................*........................",
             "...................*........................"],

            ["............................................",  # Expected genotype
             "...................*........................"]
        )

    def test_calls_homozygous_deletion(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            ["...................*........................",
             "...................*........................",
             "...................*........................",
             "...................*........................"],

            ["...................*........................",  # Expected genotype
             "...................*........................"]
        )


class TestCallingCombinationsOfVariants(AsciiWecallRunnerTest):

    def test_calls_heterozygous_snps_on_both_edges_further_than_15_bases(self):
        self.calls_variants_with_genotype(
            "ACGCCCCCATGCTCCAAAGAA",
            ["....................C",
             "....................C",
             "C....................",
             "C...................."],

            ["....................C",
             "C...................."]  # Expected genotype
        )

    def test_calls_snps_on_both_edges_when_gap_at_least5_low_coverage(self):
        self.calls_variants_with_genotype(
            "ACGCCCCCTGCAAAAAAAAAAAAAAAAAAA",
            ["G.......................      ",
             "G.......................      ",
             "G.............................",
             ".............................C",
             ".............................C",
             ".............................C"],

            ["G.............................",
             ".............................C"]  # Expected genotype
        )

    def test_calls_two_heterozygous_snps_at_same_pos(self):
        svc_driver = SVCDriver(self)
        svc_driver\
            .with_ref_sequence(
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom="1")\
            .with_read(
                "...................T........................", n_fwd=1, n_rev=1, chrom="1", sample_name="sample")\
            .with_read(
                "...................C........................", n_fwd=1, n_rev=1, chrom="1", sample_name="sample")

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect\
            .has_record("1", 19, "G", "T")\
            .with_sample("sample")\
            .has_genotype("./1")
        vcf_expect\
            .has_record("1", 19, "G", "C")\
            .with_sample("sample")\
            .has_genotype("./1")

    def test_calls_two_heterozygous_snps_at_different_pos(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            ["...................T........................",
             "...................T........................",
             "....................C.......................",
             "....................C......................."],

            ["...................T........................",  # Expected genotype
             "....................C......................."]
        )

    def test_calls_two_heterozygous_deletions_at_same_pos(self):
        svc_driver = SVCDriver(self)
        svc_driver\
            .with_ref_sequence(
                "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG", chrom="1")\
            .with_read(
                "...................**.......................", n_fwd=1, n_rev=1, chrom="1", sample_name="sample")\
            .with_read(
                "...................*........................", n_fwd=1, n_rev=1, chrom="1", sample_name="sample")

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect\
            .has_record("1", 18, "AG", "A")\
            .with_sample("sample")\
            .has_genotype("./1")
        vcf_expect\
            .has_record("1", 18, "AGT", "A")\
            .with_sample("sample")\
            .has_genotype("./1")

    def test_calls_two_heterozygous_deletions_at_different_pos(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCATG",
            ["...................**.......................",
             "...................**.......................",
             ".....................*......................",
             ".....................*......................"],

            ["...................**.......................",  # Expected genotype
             ".....................*......................"]
        )

    def test_calls_two_heterozygous_insertions_at_same_pos(self):
        svc_driver = SVCDriver(self)
        svc_driver\
            .with_ref_sequence(
                "AAAGCGTACAACCGGGTTC**GTCACAAACCCGTTACGTATATG", chrom="1")\
            .with_read(
                "...................AA.......................", n_fwd=1, n_rev=1, chrom="1", sample_name="sample")\
            .with_read(
                "...................*A.......................", n_fwd=1, n_rev=1, chrom="1", sample_name="sample")

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect\
            .has_record("1", 18, "C", "CA")\
            .with_sample("sample")\
            .has_genotype("./1")
        vcf_expect\
            .has_record("1", 18, "C", "CAA")\
            .with_sample("sample")\
            .has_genotype("./1")

    def test_calls_two_heterozygous_insertions_at_different_pos(self):
        self.calls_variants_with_genotype(
            "AAAGCGTACAACCGGGTTC**GTCACAAAC**CCGTCGTATATG",
            ["...................AA.........**............",
             "...................AA.........**............",
             "...................**.........AA............",
             "...................**.........AA............"],

            ["...................AA.........**............",  # Expected genotype
             "...................**.........AA............"]
        )


class TestMNPCalling(AsciiWecallRunnerTest):
    def calls_variants_with_coverage_20(self, ref, sequence_list, expected_variants=None):
        self.calls_variants(
            ref,
            sequence_list,
            expected_ascii_haplotypes=None,
            expected_variant_stubs=expected_variants,
            n_fwd=20,
            n_rev=20,
            config_dict={
                "allowMNPCalls": "True",
                "maxClusterDist": 20})

    def test_calls_mnp_on_adjacent_bases_single_read(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".......AG....................."],
            expected_variants={(7, "CT", "AG")}
        )

    def test_calls_mnp_with_bases_apart_by_1(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".......A.A...................."],
            expected_variants={(7, "CTG", "ATA")}
        )

    def test_calls_mnp_with_bases_apart_by_10(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA",
            [".......A....................A............"],
            expected_variants={(7, "CTGGGGGGGGGTGGGGGGGGGG", "ATGGGGGGGGGTGGGGGGGGGA")}
        )

    def test_calls_two_snps_with_bases_apart_by_21(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA",
            [".......A.....................A..........."]
        )

    def test_calls_mnp_with_three_snps_and_bases_apart_by_1(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".......A.A.A.................."],
            expected_variants={(7, "CTGGG", "ATAGA")}
        )

    def test_calls_mnp_with_snps_apart_by_1(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".T.G.G.A.A.A.T.T.T.A.G.G.G.G.T"],
            expected_variants={(1, "CGCCCCCTGGGGGGGGGGCAAAAAAAAAA", "TGGCGCATAGAGTGTGTGAAGAGAGAGAT")}
        )

    def test_calls_mnp_with_snps_apart_by_1_starting_in_the_middle(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".......A.A.A.T.T.T.A.G.G.G.G.T"],
            expected_variants={(7, "CTGGGGGGGGGGCAAAAAAAAAA", "ATAGAGTGTGTGAAGAGAGAGAT")}
        )

    def test_calls_mnp_with_adjacent_snps(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".TTGTAGAGC...................."],
            expected_variants={(1, "CGCCCCCTG", "TTGTAGAGC")}
        )

    def test_calls_mnp_with_competely_different_seq(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            ["TTTGTAGAGCCATCATCATATTTGGGCCTG"],
            expected_variants={(0, "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA", "TTTGTAGAGCCATCATCATATTTGGGCCTG")}
        )

    def test_calls_adjacent_snps_and_mnp_on_separate_reads(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".......A.T....................",
             "........G....................."],
            expected_variants={(7, "CTG", "ATT"), (8, "T", "G")}
        )

    def test_calls_overlapping_mnps_on_separate_reads(self):
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGGGGCAAAAAAAAAA",
            [".......A.T......T.............",
             "....A...G...C......T...C......"],
            {(7, "CTGGGGGGGG", "ATTGGGGGGT"), (4, "CCCCTGGGGGGGGGGCAAAA", "ACCCGGGGCGGGGGGTAAAC")}
        )

    def test_calls_mnps_when_interupted_by_deletion_within_too_long_homopolymer(self):
        # deletion not called.  If the GG homopolymer one shorter it does get called
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGAGGGGGGGGGGGGGCAAAAAAAAAA",
            [".......A.A.A.T.T.*............A.G.G.G.G.T"],
            {(7, "CTGGGGGGG", "ATAGAGTGT"), (16, "AG", "A"), (30, "CAAAAAAAAAA", "AAGAGAGAGAT")}
        )

    def test_calls_mnps_when_interupted_by_deletion_shorter_homopolymer(self):
        # deletion now called because homopolymer just short enough
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGAGGGGGGGGGGGCAAAAAAAAAA",
            [".......A.A.A.T.T.*..........A.G.G.G.G.T"],
            {(7, "CTGGGGGGG", "ATAGAGTGT"), (16, "AG", "A"), (28, "CAAAAAAAAAA", "AAGAGAGAGAT")}
        )

    def test_calls_del_when_followed_by_mnp(self):
        # same example as above (test_calls_mnps_when_interupted_by_deletion_within_too_long_homopolymer)
        # but the first MNP removed.  Deletion still not called
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGAGGGGGGGGGGGGGCAAAAAAAAAA",
            [".................*............A.G.G.G.G.T"],
            {(16, "AG", "A"), (30, "CAAAAAAAAAA", "AAGAGAGAGAT")}
        )

    def test_calls_del_when_preceded_by_mnp(self):
        # same example as above (test_calls_mnps_when_interupted_by_deletion_within_too_long_homopolymer)
        # but the second MNP has now been removed.  Deletion now called
        self.calls_variants_with_coverage_20(
            "ACGCCCCCTGGGGGGGAGGGGGGGGGGGGGCAAAAAAAAAA",
            [".......A.A.A.T.T.*......................."],
            {(7, "CTGGGGGGG", "ATAGAGTGT"), (16, "AG", "A")}
        )


class TestVariantCallingWithCustomQuality(AsciiWecallRunnerTest):
    def test_calls_two_snps_with_high_quality(self):
        self.calls_variants(
            "ACGCCCCCTGCAAAAAAAAAA",  # input
            [".T...................",
             " 9                   ",
             ".T...................",
             " 9                   ",
             "...........C.........",
             "           9         ",
             "...........C.........",
             "           9         "]
        )

    def test_should_not_call_snp_with_low_quality(self):
        self.calls_variants(
            "ACGCCCCCTGCAAAAAAAAAA",  # input
            [".T...................",
             " 1                   ",
             ".T...................",
             " 1                   ",
             "...........C.........",
             "           9         ",
             "...........C.........",
             "           9         "],

            ["...........C.........",  # expected output
             "....................."]
        )

    def test_should_not_call_snp_with_low_quality_in_the_middle(self):
        self.calls_variants(
            "ACGCCCCCTGCAAAAAAAAAA",  # input
            [".T...................",
             " 9                   ",
             ".T...................",
             " 9                   ",
             "...........C.........",
             "           1         ",
             "...........C.........",
             "           1         "],

            [".T...................",  # expected output
             "....................."]
        )

    def test_should_not_call_snp_with_low_quality_on_different_length_reads(self):
        self.calls_variants(
            "ACGCCCCCTGCAAAAAAAAAACCCCCCTTTGGGGGGGGGG",  # input
            ["     .T......................           ",
             "      1                                 ",
             "......T...........................      ",
             "      1                                 ",
             "          ...........G..................",
             "                     9                  ",
             "  ...................G.........         ",
             "                     9                  "],

            ["........................................",  # expected output
             ".....................G.................."]
        )


class TestCallingInExtremeEdgeCases(AsciiWecallRunnerTest):
    def test_should_not_call_anything_in_silly_case(self):
        self.calls_variants(
            "A**********************************************A",
            [".**********************************************.",
             ".**********************************************."]
        )
