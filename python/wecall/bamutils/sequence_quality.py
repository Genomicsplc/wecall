# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.sequence_position import MISSING_BASE
import re
from wecall.common.exceptions import weCallException


class SequenceQuality(object):
    QUALITY_MAPPING = {
        "0": 0,
        "1": 10,
        "2": 20,
        "3": 30,
        "4": 40,
        "5": 50,
        "6": 60,
        "7": 70,
        "8": 80,
        "9": 90
    }

    DEFAULT_QUALITY = 39

    def __init__(self, quality_string, quality_mapping=QUALITY_MAPPING):
        if not SequenceQuality.is_valid_qual(quality_string):
            raise weCallException(
                "Illegal character in the quality string {!r}".format(quality_string))

        self.quality_mapping = quality_mapping
        self.ascii_quality = self.parse_quality_to_ascii(quality_string)

    def parse_quality_to_ascii(self, quality_string):
        ascii_quality = ""
        for qual_char in quality_string:
            if qual_char == MISSING_BASE:
                ascii_quality += SequenceQuality.to_ascii(
                    SequenceQuality.DEFAULT_QUALITY)
            else:
                ascii_quality += SequenceQuality.to_ascii(
                    self.quality_mapping[qual_char])

        return ascii_quality

    @staticmethod
    def parse_ascii_to_quality_string(ascii_quality):
        quality_string = ""
        for ascii_char in ascii_quality:
            quality = SequenceQuality.from_ascii(ascii_char)
            if quality == SequenceQuality.DEFAULT_QUALITY:
                quality_string += MISSING_BASE
            else:
                quality_string += str(int(quality / 10))
        return quality_string

    @staticmethod
    def to_ascii(quality_score):
        return chr(quality_score + ord('!'))

    @staticmethod
    def from_ascii(char_quality):
        return ord(char_quality) - ord('!')

    @staticmethod
    def is_valid_qual(quality_string):
        return re.match(r'^[ \d]*\Z', quality_string)
