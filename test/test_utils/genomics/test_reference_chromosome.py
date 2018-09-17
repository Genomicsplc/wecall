# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from wecall.common.exceptions import EchidnaException
from wecall.genomics.reference_chromosome import ReferenceChromosome


class TestSequenceReference(TestCase):
    def test_should_raise_when_gap_in_reference(self):
        self.assertRaisesRegex(
            EchidnaException,
            "Illegal character in reference sequence \'TA AA\'",
            ReferenceChromosome,
            "TA AA"
        )

    def test_should_raise_when_unknown_character_in_reference(self):
        self.assertRaisesRegex(
            EchidnaException,
            "Illegal character in reference sequence \'TA&AA\'",
            ReferenceChromosome,
            "TA&AA"
        )

    def test_should_raise_when_dot_in_reference(self):
        self.assertRaisesRegex(
            EchidnaException,
            "Illegal character in reference sequence \'TA.AA\'",
            ReferenceChromosome,
            "TA.AA"
        )

    def test_should_allocate_pos_from_and_pos_to_based_on_reference_size(self):
        # Given
        input_ref_seq = "ACCCT"
        # When
        seq_ref = ReferenceChromosome(input_ref_seq)
        # Then
        self.assertEqual(seq_ref.pos_from, 0)
        self.assertEqual(seq_ref.pos_to, len(input_ref_seq))

    def test_should_ignore_asterixes_in_reference_sequence_in_computing_pos_to(self):
        # Given
        ref_seq = ReferenceChromosome("C*C*C")
        # Then
        self.assertEqual(ref_seq.pos_to, 3)

    def test_should_be_able_to_access_ref_char(self):
        # Given
        seq_ref = ReferenceChromosome("AC*T*G")
        # Then
        self.assertEqual(seq_ref[0], "A")
        self.assertEqual(seq_ref[1], "C")
        self.assertEqual(seq_ref[2], "T")
        self.assertEqual(seq_ref[3], "G")

    def test_should_correctly_getitem_for_offset_reference(self):
        # Given
        seq_ref = ReferenceChromosome("A*T", 10)
        # Then
        self.assertEqual(seq_ref[10], "A")
        self.assertEqual(seq_ref[11], "T")

    def test_should_get_correct_fasta_string(self):
        # Given
        seq_ref = ReferenceChromosome("AC*T*G")
        # Then
        self.assertEqual(seq_ref.fasta_string(), "ACTG")

    def test_should_get_correct_fasta_string_for_offset_reference(self):
        # Given
        seq_ref = ReferenceChromosome("AC*T*G", 5)
        # Then
        self.assertEqual(seq_ref.fasta_string(), "NNNNNACTG")
