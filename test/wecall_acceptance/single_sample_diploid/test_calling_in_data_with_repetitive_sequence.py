# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestLeftAlignment(AsciiWecallRunnerTest):
    def test_calls_deletion_aligned_to_homopolymer_start(self):
        self.calls_variants(
            "AAATGAAAAAAAAAACTGTTACGGC",
            [".............*...........",
             "............*............",
             ".........*..............."],
            [".....*...................",
             ".....*..................."]
        )

    def test_calls_deletion_aligned_to_complex_polymer_repeat_units(self):
        self.calls_variants(
            "AAATGATCGTATCGTATCGTATCGTATCGTATCGTATCGTATCGTG",
            ["..................................*****.......",
             ".................................*****........",
             "..............................*****..........."],
            [".....*****....................................",
             ".....*****...................................."]
        )

    def test_calls_deletion_correctly_after_a_snp(self):
        svc_driver = SVCDriver(self) \
            .with_ref_sequence(
                "AAATGAAAAAAAAAACTGTTACGGC", chrom="1"
        ).with_read(
                "......T......*...........", chrom="1"
        ).with_read(
                "......T.....*............", chrom="1"
        ).with_read(
                "......T..*...............", chrom="1"
        )

        expect = svc_driver.call().with_output_vcf()

        expect.record_count(2)

        expect.has_record("1", 4, "GA", "G")
        expect.has_record("1", 7, "A", "T")

    def test_calls_deletion_correctly_around_a_snp(self):
        # The SNP should mean that the indel could be left-aligned further.
        svc_driver = SVCDriver(self) \
            .with_ref_sequence(
                "AAAGAATAAAAAAAACTGTTACGGC", chrom="1"
        ).with_read(
                "......A......*...........", chrom="1"
        ).with_read(
                "......A.....*............", chrom="1"
        ).with_read(
                "......A*.................", chrom="1"
        )

        expect = svc_driver.call().with_output_vcf()

        expect.record_count(1)

        expect \
            .has_record("1", 5, "AT", "A")


class TestSingleDeletionCallingInRepetitiveSequence(AsciiWecallRunnerTest):

    def test_calls_deletion_to_left_of_A10_homopolymer(self):
        self.calls_variants(
            "AAATGAAAAAAAAAACTGTTACGGC",
            ["....*....................",
             "....*...................."]
        )

    def test_calls_deletion_to_left_of_A15_homopolymer(self):
        self.calls_variants(
            "AAATGAAAAAAAAAAAAAAATCGGC",
            ["....*....................",
             "....*...................."]
        )

    def test_calls_deletion_to_left_of_A18_homopolymer(self):
        self.calls_variants(
            "AAATGAAAAAAAAAAAAAAAAAAGC",
            ["....*....................",
             "....*...................."]
        )

    def test_calls_deletion_to_left_of_A20_homopolymer(self):
        self.calls_variants(
            "AAATGAAAAAAAAAAAAAAAAAAAACTC",
            ["....*.......................",
             "....*.......................",
             "....*......................."]
        )

    def test_calls_deletion_to_left_of_A20_homopolymer_with_right_anchor(self):
        self.calls_variants(
            "AAATGAAAAAAAAAAAAAAAAAAAACT",
            ["....*......................",
             "....*......................",
             "....*......................"]
        )


class TestSNPAndDeletionCallingInRepetitiveSequence(AsciiWecallRunnerTest):

    def test_calls_correct_overlapping_deletion_and_snp_in_non_repetitive_sequence(self):
        self.calls_variants(
            "TGTCAGGACATGGCATAACAAGATAC",
            ["......T...................",
             "......T...................",
             ".....**...................",
             ".....**..................."]
        )

    def test_calls_correct_overlapping_deletion_and_snp_in_A2_homopolymer(self):
        self.calls_variants(
            "TGTCGTAACATGGCATAACAAGATAC",
            ["......T...................",
             "......T...................",
             ".....**...................",
             ".....**..................."]
        )

    def test_calls_correct_overlapping_deletion_and_snp_in_A3_homopolymer(self):
        self.calls_variants(
            "TGTCGAAACATGGCATAACAAGATAC",
            ["......T...................",
             "......T...................",
             ".....**...................",
             ".....**..................."]
        )

    def test_calls_correct_overlapping_deletion_and_snp_in_A4_homopolymer(self):
        svc_driver = SVCDriver(self) \
            .with_ref_sequence(
                "TGTCAAAACATGGCATAACAAGATAC", chrom="1"
        ).with_read(
                "......T...................", chrom="1", n_rev=1, n_fwd=1
        ).with_read(
                ".....**...................", chrom="1", n_rev=1, n_fwd=1
        )

        expect = svc_driver.call().with_output_vcf()

        expect.record_count(2)

        expect \
            .has_record("1", 3, "CAA", "C")

        expect \
            .has_record("1", 6, "A", "T")

    def test_calls_correct_overlapping_deletion_and_snp_in_A6_homopolymer(self):
        svc_driver = SVCDriver(self) \
            .with_ref_sequence(
                "TGAAAAAACATGGCATAACAAGATAC", chrom="1"
        ).with_read(
                "......T...................", chrom="1", n_rev=1, n_fwd=1
        ).with_read(
                ".....**...................", chrom="1", n_rev=1, n_fwd=1
        )

        expect = svc_driver.call().with_output_vcf()

        expect.record_count(2)

        expect \
            .has_record("1", 1, "GAA", "G")

        expect \
            .has_record("1", 6, "A", "T")


class TestCallingOverlappingDeletions(AsciiWecallRunnerTest):

    def test_calls_overlapping_deletions_in_non_repetitive_sequence(self):
        self.calls_variants(
            "CTAGAATTCCGATACAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A3_homopolymer(self):
        self.calls_variants(
            "CTAGAAATCCGATACAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A4_homopolymer(self):
        self.calls_variants(
            "CTAGAAAACCGATACAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A5_homopolymer(self):
        self.calls_variants(
            "CTAGAAAAACGATACAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A6_homopolymer(self):
        self.calls_variants(
            "CTAGAAAAAAGATACAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A7_homopolymer(self):
        self.calls_variants(
            "CTAGAAAAAAATCTCAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A8_homopolymer(self):
        self.calls_variants(
            "CTAGAAAAAAATCTCAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A9_homopolymer(self):
        self.calls_variants(
            "CTAGAAAAAAAACTCAGATAACAAACCC",
            ["....*.......................",
             "....*.......................",
             "....**......................",
             "....**......................"],
        )

    def test_calls_overlapping_deletions_in_A10_homopolymer(self):
        svc_driver = SVCDriver(self) \
            .with_ref_sequence(
                "CTAGAAAAAAAAAACAGATAACAAACCC", chrom="1"
        ).with_read(
                "....*.......................", chrom="1", n_rev=2, n_fwd=1
        ).with_read(
                "....**......................", chrom="1", n_rev=2, n_fwd=1
        )

        expect = svc_driver.call().with_output_vcf()

        expect.record_count(2)
        expect.has_record("1", 3, "GAA", "G")
        expect.has_record("1", 3, "GA", "G")

    def test_should_be_able_to_call_distinct_homopolymer_insertions_in_a_long_homopolymer_region(self):
        self.calls_variants(
            "ATTCAGATACTTTGCCCATTTTTAAGTTGGATCATTAGATTTTTTTCCTATAGAATTG****TTTTTTTTTTTTATTTCCTGTTATTAATCCCTTGTCAGATTTTTTTTTTGCAAATATTTTTT",  # noqa
            ["        ..................................................TTTT............................................................  ",  # noqa
             " .........................................................TT**...................................................           "],  # noqa
            ["..........................................................TTTT..............................................................",  # noqa
             "..........................................................TT**.............................................................."],  # noqa
            n_fwd=7, n_rev=0
        )
