# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest


class TestLeftAlignOutOfRead(AsciiWecallRunnerTest):

    def test_should_not_call_unleftalignable_insertion(self):
        self.calls_variants(
            "CTAGAAAAAAAAAAAA*AATTAAAAAAAAAAACAGATAACAAACCC",
            ["            ....A..........................   ",
             "          ......A.......................      ",
             "            ,,,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "       ,,,,,,,,,a,,,,,,,,,,,,,,,,,,,,,        "],
            expected_variant_stubs=[]
        )

    def test_should_call_left_aligned_insertion_when_inside_read(self):
        self.calls_variants(
            "CTAGAAAAAAAAAAAA*AATTAAAAAAAAAAACAGATAACAAACCC",
            ["................A..........................   ",
             "................A.......................      ",
             ",,,,,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,,,,,,,,,a,,,,,,,,,,,,,,,,,,,,,        "],
            expected_variant_stubs=[(3, "G", "GA")]
        )

    def test_should_not_call_unleftalignable_deletion(self):
        self.calls_variants(
            "CTAGAAAAAAAAAAAAAAATTAAAAAAAAAAACAGATAACAAACCC",
            ["            ....*..........................   ",
             "          ......*.......................      ",
             "            ,,,,*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "       ,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,        "],
            expected_variant_stubs=[]
        )

    def test_should_call_left_aligned_deletion_when_inside_read(self):
        self.calls_variants(
            "CTAGAAAAAAAAAAAAAAATTAAAAAAAAAAACAGATAACAAACCC",
            ["................*..........................   ",
             "................*.......................      ",
             ",,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             ",,,,,,,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,        "],
            expected_variant_stubs=[(3, "GA", "G")]
        )

    def test_should_not_call_insertion_at_start_of_reads(self):
        self.calls_variants(
            "CTAGCTGACA*ACTATTAAAAAAAAAAAAACAGATAACAAACCCCCCC",
            ["          T..................................   ",
             "          T...............................      ",
             "          t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "          t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "],
            expected_variant_stubs=[]
        )

    def test_should_not_call_insertion_at_end_of_reads(self):
        self.calls_variants(
            "CTAGCTGACAACTATTAAAAAAAAAAAAACAGATAACAA*ACCCCCCC",
            ["      .................................T        ",
             "        ...............................T        ",
             "         ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t        ",
             "          ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t        "],
            expected_variant_stubs=[]
        )

    def test_should_not_call_deletion_at_start_of_reads(self):
        self.calls_variants(
            "CTAGCTGACATACTATTAAAAAAAAAAAAACAGATAACAAACCCCCCC",
            ["          *..................................   ",
             "          *...............................      ",
             "          *,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "          *,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "],
            expected_variant_stubs=[]
        )

    def test_should_not_call_deletion_at_end_of_reads(self):
        self.calls_variants(
            "CTAGCTGACAACTATTAAAAAAAAAAAAACAGATAACAATACCCCCCC",
            ["      .................................*        ",
             "        ...............................*        ",
             "         ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*        ",
             "          ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*        "],
            expected_variant_stubs=[]
        )

    def test_should_call_snp_at_start_of_reads(self):
        self.calls_variants(
            "CTAGCTGACAGACTATTAAAAAAAAAAAAACAGATAACAAACCCCCCC",
            ["          T..................................   ",
             "          T...............................      ",
             "          t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
             "          t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,       "],
            expected_variant_stubs=[(10, "G", "T")]
        )
