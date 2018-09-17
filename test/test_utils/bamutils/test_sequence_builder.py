# All content Copyright (C) 2018 Genomics plc
import unittest
from wecall.bamutils.sequence_builder import sequence_builder
from wecall.common.exceptions import EchidnaException
from wecall.genomics.reference_chromosome import ReferenceChromosome


class TestSequenceBuilder(unittest.TestCase):

    def test_should_raise_for_invalid_char_in_seq(self):
        with self.assertRaisesRegex(EchidnaException, "Illegal character in sequence .*'"):
            sequence_builder(ReferenceChromosome("TAAAA"), "..&..")

    def test_should_raise_for_lower_case_char_in_fwd_seq(self):
        with self.assertRaisesRegex(EchidnaException, "Illegal character in sequence .*"):
            sequence_builder(ReferenceChromosome("TAAAA"), "..c..")

    def test_should_raise_for_dot_in_reverse_seq(self):
        with self.assertRaisesRegex(EchidnaException, "Illegal character in sequence .*"):
            sequence_builder(ReferenceChromosome("TAAAA"), ",,c.,")

    def test_should_raise_when_quality_string_too_short(self):
        with self.assertRaisesRegex(EchidnaException, "Quality string has to be of the same length as reference."):
            sequence_builder(ReferenceChromosome("TAAAA"), ".....", "22  ")

    def test_should_raise_when_quality_string_too_short_due_to_insertions(self):
        with self.assertRaisesRegex(EchidnaException, "Quality string has to be of the same length as reference."):
            sequence_builder(ReferenceChromosome("TA**A"), "..TT.", "1234")

    def test_should_raise_when_quality_string_too_short_multisequence(self):
        with self.assertRaisesRegex(EchidnaException, "Quality string has to be of the same length as reference."):
            sequence_builder(
                ReferenceChromosome("TAAAA*A"),
                "...  ..",
                "12   4")

    def test_should_raise_when_quality_assigned_to_gap(self):
        with self.assertRaisesRegex(EchidnaException, "Cannot assign base quality inside a gap."):
            sequence_builder(ReferenceChromosome(
                "TAAAA*A"), "...  ..", "12  34 ")

    def test_should_raise_when_quality_string_too_long_due_to_insertions(self):
        with self.assertRaisesRegex(EchidnaException, "Quality string has to be of the same length as reference."):
            sequence_builder(ReferenceChromosome("TA**A"), "..TT.", "123")

    def test_should_build_correct_sequence_without_any_whitespace(self):
        ref = ReferenceChromosome("C*CC")
        annotated_seqs = sequence_builder(ref, ".*.T")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].pos, 0)
        self.assertEqual(reads[0].rlen, 3)
        self.assertEqual(reads[0].seq, "CCT")

    def test_should_build_correct_sequence_with_insertion_at_the_end(self):
        ref = ReferenceChromosome("CCC**")
        builders = sequence_builder(ref, "...TT")
        read_lists = [builder.build_reads(0, {}) for builder in builders]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(reads[0].pos, 0)
        self.assertEqual(reads[0].rlen, 5)
        self.assertEqual(reads[0].seq, "CCCTT")

    def test_should_interpret_leading_whitespace_to_override_pos_from(self):
        ref = ReferenceChromosome("CATG")
        annotated_seqs = sequence_builder(ref, "  .T")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].pos, 2)
        self.assertEqual(reads[0].seq, "TT")

    def test_should_interpret_leading_whitespace_to_override_pos_from_when_ref_has_deletion(self):
        ref = ReferenceChromosome("C*TG")
        annotated_seqs = sequence_builder(ref, "  C.")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].pos, 1)
        self.assertEqual(reads[0].seq, "CG")

    def test_should_interpret_trailing_whitespace_to_override_pos_to(self):
        ref = ReferenceChromosome("CATG")
        annotated_seqs = sequence_builder(ref, ".C  ")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].pos, 0)
        self.assertEqual(reads[0].seq, "CC")

    def test_should_interpret_trailing_whitespace_to_override_pos_to_when_seq_has_insertion(self):
        ref = ReferenceChromosome("C*GA")
        annotated_seqs = sequence_builder(ref, ".C  ")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(reads[0].rlen, 2)

    def test_should_interpret_trailing_whitespace_to_override_positions_for_complex_ref_and_seq(self):
        ref = ReferenceChromosome("ACCC*G*A")
        annotated_seqs = sequence_builder(ref, ".**.C   ")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(reads[0].pos, 0)
        self.assertEqual(reads[0].rlen, 3)

    def test_should_translate_reverse_seq_into_correct_annotations(self):
        ref = ReferenceChromosome("CCTG")
        annotated_seq = sequence_builder(ref, ",,c,")[0]
        self.assertEqual(annotated_seq.n_fwd, 0)
        self.assertEqual(annotated_seq.n_rev, 1)

    def test_should_translate_reverse_seq_into_correct_sequence(self):
        ref = ReferenceChromosome("AAACCTG*TAA")
        builders = sequence_builder(ref, "   ,,c,*,  ")
        read_lists = [builder.build_reads(0, {}) for builder in builders]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(reads[0].seq, 'CCCGT')


class TestQualityBuilding(unittest.TestCase):

    def setUp(self):
        self.ascii_codes = {
            "0": "!",
            "1": "+",
            "2": "5",
            "3": "?",
            "4": "I",
            "5": "S",
            "6": "]",
            "7": "g",
            "8": "q",
            "9": "{"}
        self.default_qual = "H"

    def test_should_build_with_default_quality_for_None(self):
        ref = ReferenceChromosome("AAAAA")
        annotated_seqs = sequence_builder(ref, ".....")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(reads[0].qual, self.default_qual * 5)

    def test_should_build_with_custom_quality_with_del(self):
        ref = ReferenceChromosome("AAAAA")
        annotated_seqs = sequence_builder(ref, "..*..", quality_string="31 00")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))
        self.assertEqual(
            reads[0].qual,
            self.ascii_codes["3"] +
            self.ascii_codes["1"] +
            self.ascii_codes["0"] *
            2)

    def test_should_build_with_custom_quality_with_ins(self):
        ref = ReferenceChromosome("AA**A")
        annotated_seqs = sequence_builder(ref, "..CC.", quality_string="31220")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(
            reads[0].qual,
            self.ascii_codes["3"] +
            self.ascii_codes["1"] +
            self.ascii_codes["2"] *
            2 +
            self.ascii_codes["0"])

    def test_should_build_with_custom_quality_and_sequence_shorter_than_reference(self):
        ref = ReferenceChromosome("AAAAAAAAAAAA")
        builders = sequence_builder(
            ref, "  ..*..     ", quality_string="  31 0      ")
        read_lists = [builder.build_reads(0, {}) for builder in builders]
        reads = [read for read_list in read_lists for read in read_list]
        self.assertEqual(
            reads[0].qual,
            self.ascii_codes["3"] +
            self.ascii_codes["1"] +
            self.ascii_codes["0"] +
            self.default_qual)

    def test_should_raise_when_assigning_qual_to_deletion(self):
        with self.assertRaisesRegex(EchidnaException, "Cannot assign base quality to a deleted base."):
            sequence_builder(ReferenceChromosome("AAAA"), ".*..", " 1  ")


class TestMultiSequenceLines(unittest.TestCase):
    def test_should_build_two_seqs_defined_on_single_line(self):
        ref = ReferenceChromosome("AAACCTGTAA")
        annotated_seqs = sequence_builder(ref, " ...  .C. ")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))
        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].pos, 1)
        self.assertEqual(reads[0].seq, "AAC")
        self.assertEqual(reads[1].pos, 6)
        self.assertEqual(reads[1].seq, "GCA")

    def test_should_build_two_complex_seqs_defined_on_single_line(self):
        ref = ReferenceChromosome("AA*CC*TGTAAGG")
        annotated_seqs = sequence_builder(ref, " .G.  ,c,*,  ")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))
        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].pos, 1)
        self.assertEqual(reads[0].seq, "AGC")
        self.assertEqual(reads[1].pos, 4)
        self.assertEqual(reads[1].seq, "TCTA")

    def test_should_build_correct_qualities_for_two_complex_seqs_defined_on_single_line(self):
        ref = ReferenceChromosome("AA*CC*TGTAAGG")
        annotated_seqs = sequence_builder(
            ref, " .G.  ,c,*,  ", "  2   1   0  ")
        read_lists = [builder.build_reads(0, {}) for builder in annotated_seqs]
        reads = [read for read_list in read_lists for read in read_list]
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))
        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].seq, "AGC")
        self.assertEqual(reads[0].pos, 1)
        self.assertEqual(reads[0].qual, 'H5H')
        self.assertEqual(reads[1].seq, "TCTA")
        self.assertEqual(reads[1].pos, 4)
        self.assertEqual(reads[1].qual, '+HH!')
