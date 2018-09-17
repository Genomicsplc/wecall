# All content Copyright (C) 2018 Genomics plc
from wecall.common.exceptions import EchidnaException
import unittest
from wecall.bamutils.sequence_position import SequencePosition


class TestSequencePosition(unittest.TestCase):
    def test_should_recognise_gap(self):
        seq_pos = SequencePosition("C", " ", " ")
        self.assertTrue(seq_pos.is_gap)

    def test_should_recognise_sequence(self):
        seq_pos = SequencePosition("A", "C", " ")
        self.assertFalse(seq_pos.is_gap)

    def test_should_fail_at_quality_inside_gap(self):
        self.assertRaisesRegex(
            EchidnaException,
            "Cannot assign base quality inside a gap.",
            SequencePosition, "A", " ", "2"
        )

    def test_should_fail_at_missing_ref_char(self):
        self.assertRaisesRegex(
            EchidnaException,
            "Missing reference character.",
            SequencePosition, " ", "A", " "
        )

    def test_should_fail_at_empty_ref_char(self):
        self.assertRaisesRegex(
            EchidnaException,
            "All characters at sequence position has to be of length 1.",
            SequencePosition, "", "C", "2"
        )

    def test_should_fail_at_empty_seq_char(self):
        self.assertRaisesRegex(
            EchidnaException,
            "All characters at sequence position has to be of length 1.",
            SequencePosition, "A", "", "2"
        )

    def test_should_fail_at_empty_qual_char(self):
        self.assertRaisesRegex(
            EchidnaException,
            "All characters at sequence position has to be of length 1.",
            SequencePosition, "A", "C", ""
        )

    def test_should_fail_at_too_log_seq_char(self):
        self.assertRaisesRegex(
            EchidnaException,
            "All characters at sequence position has to be of length 1.",
            SequencePosition, "A", "CT", " "
        )

    def test_should_fail_at_qual_assignment_to_deleted_base(self):
        self.assertRaisesRegex(
            EchidnaException,
            "Cannot assign base quality to a deleted base.",
            SequencePosition, "A", "*", "2"
        )

    def test_should_increase_ref_position_for_matching_base(self):
        seq_pos = SequencePosition("A", "C", " ")
        ref_pos = seq_pos.update_ref_pos(2)
        self.assertEqual(ref_pos, 3)

    def test_should_increase_ref_position_for_deletion(self):
        seq_pos = SequencePosition("A", "*", " ")
        ref_pos = seq_pos.update_ref_pos(2)
        self.assertEqual(ref_pos, 3)

    def test_should_not_increase_ref_position_for_insertion(self):
        seq_pos = SequencePosition("*", "C", " ")
        ref_pos = seq_pos.update_ref_pos(2)
        self.assertEqual(ref_pos, 2)
