# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest


class TestCallingWithNsInReference(AsciiWecallRunnerTest):
    def test_should_not_call_deletion_of_gap_character(self):
        self.calls_variants(
            "TTTTTTTTTTTNTTTTTTTTTTTTTTTT",
            ["...........*................",
             "...........*................",
             "...........*................",
             "...........*................"],

            ["............................",
             "............................"],  # Expected genotype
        )

    def test_should_not_deleletion_if_ref_is_N_after_aligning(self):
        self.calls_variants(
            "TTTTTTTTTNAATCGTAATTTGACACAT",
            ["...........*................",
             "...........*................",
             "...........*................",
             "...........*................"],

            ["............................",
             "............................"],  # Expected genotype
        )

    def test_should_not_insertion_if_ref_is_N_after_aligning(self):
        self.calls_variants(
            "TTTTTTTTNTG**ATCGTAATTTGACACAT",
            ["...........TG.................",
             "...........TG.................",
             "...........TG.................",
             "...........TG................."],

            ["...........**.................",
             "...........**................."],  # Expected genotype
        )

    def test_should_not_snp_if_ref_is_N(self):
        self.calls_variants(
            "CTCNNNNNNNNNNTTTTTTTTTTTTTTT",
            ["....A.......................",
             "....A.......................",
             "....A.......................",
             "....A......................."],

            ["............................",
             "............................"],  # Expected genotype
        )


class TestCallingWithNsInReads(AsciiWecallRunnerTest):
    def test_should_not_insertion_is_N(self):
        self.calls_variants(
            "CTCT**TTTTTTTTTTCCCCCCCCCCCC",
            ["....NN......................",
             "....NN......................",
             "....NN......................",
             "....NN......................"],

            ["....**......................",
             "....**......................"],  # Expected genotype
        )

    def test_should_not_snp_if_alt_is_N(self):
        self.calls_variants(
            "CTCTTTTTTTTTTTTTTTTTTTTTTTTT",
            ["....N.......................",
             "....N.......................",
             "....N.......................",
             "....N.......................",
             "............................"],

            ["............................",
             "............................"],  # Expected genotype
        )
