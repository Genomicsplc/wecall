# All content Copyright (C) 2018 Genomics plc
from wecall.common.exceptions import EchidnaException
import unittest
from wecall.bamutils.sample_bank import SampleBank
from wecall.genomics.variant import Variant


class TestSampleBank(unittest.TestCase):
    def setUp(self):
        self.default_char = "H"

    def test_should_add_sequences_with_same_reference(self):
        sample_bank = SampleBank("AAATTTTGGGGG")
        sample_bank.add_sample_name("SAMPLE1")
        sample_bank.add_sample_name("SAMPLE2")

        self.assertEqual(
            sample_bank["SAMPLE1"].reference.ref_seq,
            sample_bank["SAMPLE2"].reference.ref_seq)

    def test_should_return_all_variants(self):
        sample_bank = SampleBank("AAATTTTGGGAG")
        sample_bank.add_sample_name("SAMPLE1")
        sample_bank.add_sample_name("SAMPLE2")

        sample_bank["SAMPLE1"].add_sequence(".....G......")
        sample_bank["SAMPLE2"].add_sequence("..........*.")

        exp_variant1 = Variant(sample_bank.reference.chrom, 5, "T", "G")
        exp_variant2 = Variant(sample_bank.reference.chrom, 9, "GA", "G")
        self.assertEqual(sample_bank["SAMPLE1"].variants, {exp_variant1})
        self.assertEqual(sample_bank["SAMPLE2"].variants, {exp_variant2})
        self.assertEqual(sample_bank.variants, {exp_variant1, exp_variant2})

    def test_should_raise_for_invalid_ref_string(self):
        self.assertRaisesRegex(
            EchidnaException,
            "Illegal character in reference sequence .*",
            SampleBank,
            "..ATTTTGGGAG"
        )

    def test_should_raise_when_adding_existing_sample(self):
        sample_bank = SampleBank("AAA")
        sample_name = "SAMPLE1"
        sample_bank.add_sample_name(sample_name)

        self.assertRaisesRegex(
            EchidnaException,
            "Sample SAMPLE1 already exists in the SampleBank.",
            sample_bank.add_sample_with_seqs_and_quals,
            sample_name,
            []
        )

    def test_should_add_sequence_with_quality(self):
        sample_bank = SampleBank("AAA")
        sample_name = "SAMPLE1"
        sample_bank.add_sample_name(sample_name)
        sample_bank[sample_name].add_sequence("...", quality_string="007")
        read_lists = [builder.build_reads(0, {})
                      for builder in sample_bank["SAMPLE1"]]
        reads = [read for read_list in read_lists for read in read_list]

        # ascii: "0": "!", "1": "+", "2": "5", "3": "?", "4": "I", "5": "S",
        # "6": "]", "7": "g", "8": "q", "9": "{"
        self.assertEqual(reads[0].qual, "!!g")

    def test_should_add_short_sequence_and_quality_list(self):
        sample_bank = SampleBank("AAA")
        sample_bank.add_sample_with_seqs_and_quals("SAMPLE1", ["...", "007"])
        read_lists = [builder.build_reads(0, {})
                      for builder in sample_bank["SAMPLE1"]]
        reads = [read for read_list in read_lists for read in read_list]

        self.assertEqual(reads[0].qual, "!!g")

    def test_should_add_seq_and_quals_list_with_deletion(self):
        sample_bank = SampleBank("AAA")
        sample_bank.add_sample_with_seqs_and_quals("SAMPLE1", [".*C", "1 3"])
        read_lists = [builder.build_reads(0, {})
                      for builder in sample_bank["SAMPLE1"]]
        reads = [read for read_list in read_lists for read in read_list]

        self.assertEqual(reads[0].qual, "+?")

    def test_should_add_two_sequence_list(self):
        sample_bank = SampleBank("AAA")
        sample_bank.add_sample_with_seqs_and_quals("SAMPLE1", ["...", " .."])
        read_lists = [builder.build_reads(0, {})
                      for builder in sample_bank["SAMPLE1"]]
        reads = [read for read_list in read_lists for read in read_list]
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))

        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].qual, "HHH")
        self.assertEqual(reads[1].qual, "HH")

    def test_should_add_two_seqs_with_one_qual_string(self):
        sample_bank = SampleBank("AAA")
        sample_bank.add_sample_with_seqs_and_quals(
            "SAMPLE1", ["...", "007", " .."])
        read_lists = [builder.build_reads(0, {})
                      for builder in sample_bank["SAMPLE1"]]
        reads = [read for read_list in read_lists for read in read_list]
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))

        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].seq, "AAA")
        self.assertEqual(reads[0].qual, "!!g")
        self.assertEqual(reads[1].seq, "AA")
        self.assertEqual(reads[1].qual, self.default_char * 2)

    def test_should_add_complex_seq_and_quals_list(self):
        sample_bank = SampleBank("AAA")
        sample_bank.add_sample_with_seqs_and_quals(
            "SAMPLE1", ["...", "007", " ..", ".*C", "1 3"])
        read_lists = [builder.build_reads(0, {})
                      for builder in sample_bank["SAMPLE1"]]
        reads = [read for read_list in read_lists for read in read_list]
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))

        self.assertEqual(len(reads), 3)
        self.assertEqual(reads[0].qual, "!!g")
        self.assertEqual(reads[1].qual, "+?")
        self.assertEqual(reads[2].qual, self.default_char * 2)

    def test_should_add_seq_and_quals_list_with_fwd_and_rev_reads(self):
        sample_bank = SampleBank("AAA")
        sample_bank.add_sample_with_seqs_and_quals(
            "SAMPLE1", ["...", "007"], n_fwd=1, n_rev=2)

        self.assertEqual(len(sample_bank["SAMPLE1"]), 1)
        self.assertEqual(sample_bank["SAMPLE1"][0].n_fwd, 1)
        self.assertEqual(sample_bank["SAMPLE1"][0].n_rev, 2)

    def test_should_raise_when_multiple_quality_strings_specified_per_sequence(self):
        sample_bank = SampleBank("AAA")

        self.assertRaisesRegex(
            EchidnaException,
            "Illegal character in sequence \'008\'",
            sample_bank.add_sample_with_seqs_and_quals,
            "SAMPLE1",
            ["...", "007", "008"]
        )

    def test_should_place_variants_at_custom_position(self):
        sample_bank = SampleBank("AAATTTTGGGAG", 100)
        sample_bank.add_sample_name("SAMPLE1")
        sample_bank.add_sample_name("SAMPLE2")

        sample_bank["SAMPLE1"].add_sequence(".....G......")
        sample_bank["SAMPLE2"].add_sequence("..........*.")

        exp_variant1 = Variant(sample_bank.reference.chrom, 105, "T", "G")
        exp_variant2 = Variant(sample_bank.reference.chrom, 109, "GA", "G")
        self.assertEqual(sample_bank["SAMPLE1"].variants, {exp_variant1})
        self.assertEqual(sample_bank["SAMPLE2"].variants, {exp_variant2})
        self.assertEqual(sample_bank.variants, {exp_variant1, exp_variant2})
