# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest


class TestAlignment(AsciiWecallRunnerTest):
    def test_should_cope_with_indel_from_two_location(self):
        self.calls_variants_with_genotype(
            "CCCCCCCTAAAAAAAAAAAAAAAAAAAAAATCCCCC",
            ["............*.......*...............",
             "............*.......*...............",
             "............*.......*...............",
             "............*.......*...............",
             "............*.......*...............",
             "............*.......*...............", ],
            ["........**..........................",
             "........**.........................."]
        )

    def test_should_cope_with_indel_with_read_errors_at_different_places_which_each_block_the_left_alignment(self):
        self.calls_variants_with_genotype(
            "CCCCCCCTAAAAAAAAAAAAAAAAAAAAAATCCCCC",
            ["...........T.......**...............",
             "...........C.......**...............",
             "...........G.......**...............",
             "............T......**...............",
             "............G......**...............",
             "............C......**...............",
             ".............T.....**...............",
             ".............C.....**...............",
             ".............G.....**...............",
             "..............T....**...............",
             "..............C....**...............",
             "..............G....**...............",
             "...............T...**...............",
             "...............C...**...............",
             "...............G...**...............",
             "................G..**...............", ],
            ["........**..........................",
             "........**.........................."]
        )

    def test_should_left_align_in_complex_repeat_pattern(self):
        self.calls_variants_with_genotype(
            "C***ATGATGATGATGATGAT*ATATA*AAA*AAC",
            [".***.................G.....T...A...",
             ".***.................G.....T...A...",
             ".***.................G.....T...A...",
             ".***.................G.....T...A...", ],
            [".ATG.................*.....*...*...",
             ".ATG.................*.....*...*..."]
        )

    def test_should_represent_neighbouring_indels_as_snp(self):
        self.calls_variants_with_genotype(
            "ATCTATGAT***GATGATGATGATATATAATGA",
            [".........TAT***..................",
             ".........TAT***..................",
             ".........TAT***..................",
             ".........TAT***..................", ],
            [".........***T....................",
             ".........***T....................", ]
        )
