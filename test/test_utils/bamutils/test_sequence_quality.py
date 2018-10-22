# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from wecall.bamutils.sequence_quality import SequenceQuality
from wecall.common.exceptions import weCallException


class TestGetVariantsFromSequence(TestCase):

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

    def test_should_raise_on_non_numerical_character(self):
        self.assertRaisesRegex(
            weCallException,
            "Illegal character in the quality string .*",
            SequenceQuality,
            "  A "
        )

    def test_should_create_default_string_for_empty_input(self):
        qual = SequenceQuality("     ")
        self.assertEqual(qual.ascii_quality, self.default_qual * 5)

    def test_should_correctly_assign_to_ascii(self):
        self.assertEqual(SequenceQuality.to_ascii(0), '!')
        self.assertEqual(SequenceQuality.to_ascii(20), '5')
        self.assertEqual(SequenceQuality.to_ascii(90), '{')

    def test_should_map_a_fully_qualified_string(self):
        qual = SequenceQuality("0123456789")
        for i in range(10):
            self.assertEqual(qual.ascii_quality[i], SequenceQuality.to_ascii(
                SequenceQuality.QUALITY_MAPPING[str(i)]))

    def test_should_correctly_map_gaps_to_default(self):
        qual = SequenceQuality(" 6 1 ")
        expected_qualities = self.default_qual + self.ascii_codes["6"] + self.default_qual \
            + self.ascii_codes["1"] + self.default_qual
        self.assertEqual(qual.ascii_quality, expected_qualities)


class TestQualityToAndFromAscii(TestCase):
    def test_should_convert_numeral(self):
        self.assertEqual(SequenceQuality.to_ascii(0), '!')

    def test_to_and_from_ascii_should_be_one_to_one_up_to_256(self):
        for i in range(0, 256 - ord('!')):
            self.assertEqual(
                i, SequenceQuality.from_ascii(
                    SequenceQuality.to_ascii(i)))


class TestUnparsingAsciiQualityString(TestCase):
    def test_should_be_able_to_convert_from_ascii_quality_string(self):

        ascii_quality_string = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"  # noqa
        quality_string = "000000000011111111112222222222333333333" + " " + "444444444455555555556666666666777777777788888888889999"  # noqa

        self.assertEqual(SequenceQuality.parse_ascii_to_quality_string(
            ascii_quality_string), quality_string)
