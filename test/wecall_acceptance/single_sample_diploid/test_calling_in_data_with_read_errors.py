# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest


class TestCallingFromDataWithReadErrors(AsciiWecallRunnerTest):

    def test_calls_heterozygous_snp_on_reads_with_scattered_base_calling_errors(self):
        self.calls_variants_with_genotype(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["       ..............C.........G...       ",
             "  ,,,c,,,t,,,,,,,,,,,c,,,,,,,,,           ",
             "    ,,,,,,,,,,,,,,,,,c,,,,,,,,,,,,,,,     ",
             "                 ....C.......G....C.......",
             "       ..........T...C.............       ",
             "  ,,t,,,,,,,,,,,,,,,,c,,,,,,,,,           "],

            [".....................C....................",  # expected output
             ".....................C...................."]
        )

    def test_calls_homozygous_snp_on_reads_with_base_calling_errors_at_every_position(self):
        self.calls_variants_with_genotype(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["       .C....G....G..C.T.G.....G...       ",
             "  ,t,c,,,t,,,,,,,,,g,c,,c,,,,,,           ",
             "    g,,,,,,t,,t,,,,,gcg,,a,,,,,,a,,g,     ",
             "                 ....C....T..G...GC..A.A.A",
             "       G..T.T....T...C.....G......G       ",
             "  c,t,t,,,,,,,,tt,,,,c,,,,,,g,g           "],

            [".....................C....................",  # expected output
             ".....................C...................."]
        )

    @expectedFailure
    def test_calls_heterozygous_snp_on_reads_with_base_calling_errors_at_every_position(self):
        self.calls_variants_with_genotype(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["       .C....G....G....T.G.....G...       ",
             "  ,t,c,,,t,,,,,,,,,g,c,,c,,,,,,           ",
             "    g,,,,,,t,,t,,,,,g,g,,a,,,,,,a,,g,     ",
             "                 ....C....T..G...GC..A.A.A",
             "       G..T.T....T...C.....G......G       ",
             "  c,t,t,,,,,,,,tt,,,,c,,,,,,g,g           "],

            ["..........................................",  # expected output
             ".....................C...................."]
        )

    def test_calls_del_on_noisy_reads(self):
        self.calls_variants(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["       ..............*.............       ",
             "  ,,,,,,,t,,,,,,,,,,,*,,,,,,,,,           ",
             "    ,*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,     ",
             "                 ............G....C.......",
             "       ..........T...*.............       ",
             "  ,,,,,,,,,,*,,,,,,,,*,,,,,,,,,           "],

            [".....................*....................",  # expected output
             ".....................*...................."]
        )


class TestCallingFromSimulationOfRealData(AsciiWecallRunnerTest):
    def test_calls_nonconflicting_variants(self):
        """
        This is a real world example, githash 4cce80b721ffc564b21682b07cbd6d4924045112
        calls a conflicting combination of a het MNP and hom SNP.  Do no regress
        back test.
        """
        self.calls_variants(
            # 012345678
            "GACCATCCCGGCTAAAACGGTGAAACCCAGTCTCTACTAAAAATACAAAA",
            [",,,,,,,, ......C............C.....................",
             "...      ,,,,,,c,,,,,,,,,,,,c,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,               .....C.....................",
             ",,,,,,,,t,,,,,,c,,,,,        .....................",
             "........T......C..A.........    ..................",
             "........T......C..A.........C.......... ..........",
             ",,,,,,,,t,,,,,,c,,,,,,,,,,,,c,,,,,,,,,,  .........",
             ",,,,,,,,t,,,,,,c,,a,,,,,,,,,c,,,,,,,,,,,  ,,,,,,,,",
             ",,,,,,,,t,,,,,,c,,a,,,,,,,,,c,,,,,,,,,,,,,,    ,,,",
             "........T......C..A.........C..................   ",
             ",,,,c,,,t,,,,,,c,,,,,,,,,,,,c,,,,,,, .............",
             ",,,,,,,,t,,,,,,c,,,,,,,,,,,,c,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,t,,,,,,c,,,,,,,,,,,,c,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,t,,,,,,c,,,,,,,,,,,,c,,,,,,,,,,,,,,,,,,,,,",
             "........T......C..A.........C.....................",
             "..................................................",
             ",,,,,,,,t,,,,,,c,,a,,,,,,,,,c,,,,,,,,,,,,,,,,,,,,,",
             "........T......C..A.........C.....................",
             "........T......C..A.........C.....................",
             "........T......C..A.........C.....................",
             "........T......C..A.........C.....................",
             ",,,,,,,,t,,,,,,c,,a,,,,,,,,,c,,,,,,,,,,,,,,,,,,,,,",
             "........T......C............C.....................",
             "........T......C..A.........C.....................",
             ",,,,,,,,t,,,,,,c,,a,,,,,,,,,c,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,t,,,,,ccc,,,,,,,,,,,c,,,,,,,,,,,,,,,,a,,,,"],

            expected_variant_stubs=[
                (8, "CGGCTAAAACGGTGAAACCCA", "TGGCTAACACAGTGAAACCCC"),
                (8, "CGGCTAAAACGGTGAAACCCA", "TGGCTAACACGGTGAAACCCC"),
            ], config_dict={"allowMNPCalls": "True"}
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
