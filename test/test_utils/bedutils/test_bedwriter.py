# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from wecall.bedutils.bedrecord import BEDRecord
from wecall.bedutils.bedwriter import bed_line_from_chrom_interval, BEDWriter
from wecall.utils.interval import ChromInterval


class MockStream(object):
    def __init__(self):
        self.lines = []

    def write(self, line):
        self.lines.append(line)


class TestBedLineFromChromInterval(TestCase):
    def test_should_write_tab_delimited_region(self):
        region = ChromInterval("20", 1, 2)

        output_line = bed_line_from_chrom_interval(region)
        self.assertEqual(output_line, "20\t1\t2")


class TestBEDWriterWritesChromIntervals(TestCase):
    def test_should_write_line_to_stream(self):
        output_stream = MockStream()

        writer = BEDWriter(output_stream)
        writer.write_chrom_interval(ChromInterval("20", 1, 2))

        # Then
        self.assertEqual(output_stream.lines, ["20\t1\t2\n"])

    def test_should_write_bed_record_to_stream(self):
        output_stream = MockStream()

        writer = BEDWriter(output_stream)
        writer.write_bed_record(BEDRecord('1', 1, 2, None, 5, 'd', 'bah'))

        self.assertEqual(output_stream.lines, ['1\t1\t2\t.\t5\td\n'])
