# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from wecall.bamutils.sequence import Sequence
from wecall.common.exceptions import EchidnaException
from wecall.genomics.reference_chromosome import ReferenceChromosome
from wecall.genomics.variant import Variant


class TestGetVariantsFromSequence(TestCase):

    def test_should_raise_when_whitespace_in_seq(self):
        ref = ReferenceChromosome("AAAA")
        self.assertRaisesRegex(
            EchidnaException,
            "Illegal character in sequence .*",
            Sequence,
            ref,
            " *.."
        )

    def test_should_raise_when_illegal_char_in_seq(self):
        ref = ReferenceChromosome("AAAA")
        self.assertRaisesRegex(
            EchidnaException,
            "Illegal character in sequence \'..*F\'",
            Sequence,
            ref,
            "..*F"
        )

    def test_should_raise_if_ref_and_seq_have_different_length(self):
        ref = ReferenceChromosome("AAA")
        self.assertRaises(EchidnaException, Sequence, ref, "..")

    def test_finds_snp(self):
        ref = ReferenceChromosome("AAAAAAAAAAAAA")
        seq = Sequence(ref, ".C...........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 1, "A", "C")})

    def test_should_find_multiple_snps(self):
        ref = ReferenceChromosome("AAAAAAAAAAAAA")
        seq = Sequence(ref, ".C.........T.")
        self.assertEqual(
            seq.variants, {
                Variant(
                    ref.chrom, 1, "A", "C"), Variant(
                    ref.chrom, 11, "A", "T")})

    def test_finds_snp_after_asterix(self):
        ref = ReferenceChromosome("T*CATAAAAAAAA")
        seq = Sequence(ref, ".*.C.........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 2, "A", "C")})

    def test_finds_snp_at_the_ref_start(self):
        ref = ReferenceChromosome("CATAAAAAAAA")
        seq = Sequence(ref, "T..........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 0, "C", "T")})

    def test_finds_snp_at_the_ref_end(self):
        ref = ReferenceChromosome("CATAAAAAAAT")
        seq = Sequence(ref, "..........C")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 10, "T", "C")})

    def test_find_single_base_deletion(self):
        ref = ReferenceChromosome("TTAAAAAAAAAAT")
        seq = Sequence(ref, "..*..........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 1, "TA", "T")})

    def test_find_multiple_single_base_deletion(self):
        ref = ReferenceChromosome("TTAAAAAGAAAAT")
        seq = Sequence(ref, "..*.....*....")
        self.assertEqual(
            seq.variants, {
                Variant(
                    ref.chrom, 1, "TA", "T"), Variant(
                    ref.chrom, 7, "GA", "G")})

    def test_find_multi_base_deletion(self):
        ref = ReferenceChromosome("TTAGCAAAAAAAT")
        seq = Sequence(ref, "..***........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 1, "TAGC", "T")})

    def test_find_multi_base_deletion_with_deletion_in_reference(self):
        ref = ReferenceChromosome("TTA*AAAAAAAAT")
        seq = Sequence(ref, "..**.........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 1, "TA", "T")})

    def test_should_not_find_deletion_on_left_edge(self):
        ref = ReferenceChromosome("TAGCAAAAAAAT")
        seq = Sequence(ref, "*...........")
        print((seq.variants))
        self.assertEqual(len(seq.variants), 0)

    def test_should_not_find_deletion_on_right_edge(self):
        ref = ReferenceChromosome("TTAGCAAAAAACT")
        seq = Sequence(ref, "............*")
        self.assertEqual(len(seq.variants), 0)

    def test_should_not_find_long_deletion_on_left_edge(self):
        ref = ReferenceChromosome("TAGCAAAAAAAT")
        seq = Sequence(ref, "***.........")
        self.assertEqual(len(seq.variants), 0)

    def test_should_not_find_long_deletion_on_right_edge(self):
        ref = ReferenceChromosome("TTAGCAAAAAACT")
        seq = Sequence(ref, "..........***")
        self.assertEqual(len(seq.variants), 0)

    def test_find_adjacent_snp_and_deletion(self):
        ref = ReferenceChromosome("TTAAAAAAAAAT")
        seq = Sequence(ref, ".G*.........")
        self.assertEqual(
            seq.variants, {
                Variant(
                    ref.chrom, 1, "T", "G"), Variant(
                    ref.chrom, 1, "TA", "T")})

    def test_find_adjacent_deletion_and_snp(self):
        ref = ReferenceChromosome("TCATAAAAAAAT")
        seq = Sequence(ref, ".*G.........")
        self.assertEqual(
            seq.variants, {
                Variant(
                    ref.chrom, 0, "TC", "T"), Variant(
                    ref.chrom, 2, "A", "G")})

    def test_should_find_single_base_insertion(self):
        ref = ReferenceChromosome("CT*AAAAAAAAAT")
        seq = Sequence(ref, "..G..........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 1, "T", "TG")})

    def test_should_find_multi_base_insertion(self):
        ref = ReferenceChromosome("CT**AAAAAAAAT")
        seq = Sequence(ref, "..GC.........")
        self.assertEqual(seq.variants, {Variant(ref.chrom, 1, "T", "TGC")})

    def test_find_adjacent_insertion_and_snp(self):
        ref = ReferenceChromosome("T*ATAAAAAAAT")
        seq = Sequence(ref, ".CG.........")
        self.assertEqual(
            seq.variants, {
                Variant(
                    ref.chrom, 0, "T", "TC"), Variant(
                    ref.chrom, 1, "A", "G")})

    def test_find_adjacent_snp_and_insertion(self):
        ref = ReferenceChromosome("TA*AAAAAAAT")
        seq = Sequence(ref, ".GC........")
        self.assertEqual(
            seq.variants, {
                Variant(
                    ref.chrom, 1, "A", "G"), Variant(
                    ref.chrom, 1, "A", "AC")})

    def test_find_multiple_variants(self):
        ref = ReferenceChromosome("TA*AAAGCTAACT")
        seq = Sequence(ref, ".GC...T...**.")
        self.assertEqual(seq.variants, {
            Variant(ref.chrom, 1, "A", "G"),
            Variant(ref.chrom, 1, "A", "AC"),
            Variant(ref.chrom, 5, "G", "T"),
            Variant(ref.chrom, 8, "AAC", "A")
        })

    def test_raise_at_dot_overlapping_asterix(self):
        ref = ReferenceChromosome("TA*AAAAAAAT")
        self.assertRaisesRegex(
            EchidnaException,
            "Invalid sequence at ref position 1.",
            Sequence,
            ref,
            "...........")


class TestGetCigarFromSequence(TestCase):
    def test_should_get_empty_cigar(self):
        ref = ReferenceChromosome("")
        seq = Sequence(ref, "")
        self.assertEqual(str(seq.cigar), "")

    def test_should_get_correct_cigar_for_dots(self):
        ref = ReferenceChromosome("CCAA")
        seq = Sequence(ref, "....")
        self.assertEqual(str(seq.cigar), "4M")

    def test_should_get_correct_cigar_for_snp(self):
        ref = ReferenceChromosome("TTT")
        seq = Sequence(ref, ".G.")
        self.assertEqual(str(seq.cigar), "3M")

    def test_should_get_correct_cigar_for_padding_and_deletion(self):
        ref = ReferenceChromosome("TTT")
        seq = Sequence(ref, ".*.")
        self.assertEqual(str(seq.cigar), "1M1D1M")

    def test_should_get_correct_cigar_for_padding_and_insertion(self):
        ref = ReferenceChromosome("T*T")
        seq = Sequence(ref, ".C.")
        self.assertEqual(str(seq.cigar), "1M1I1M")

    def test_should_get_correct_cigar_for_multiple_events(self):
        ref = ReferenceChromosome("CCC***AAATTT")
        seq = Sequence(ref, "A*.T*T...**C")
        self.assertEqual(str(seq.cigar), "1M1D1M2I3M2D1M")
