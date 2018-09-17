# All content Copyright (C) 2018 Genomics plc
import os
import re
import unittest

from wecall.genomics.variant import Variant
from wecall.vcfutils.genotype_call import GenotypeCall
from wecall.vcfutils.parser import VCFReader, VCFReaderContextManager, decode_VCF_string, \
    parse_VCF_comma_separated_pair_value
from wecall.vcfutils.schema import Schema
from wecall.vcfutils.writer import VCFWriterContextManager
from wecall_test_drivers.base_test import BaseTest


class ParserTest(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.data_dir = os.path.join(os.path.dirname(__file__), "example_data")

    def variant_is_equal(self, var1, var2):
        self.assertEqual(var1.chrom, var2[0])
        self.assertEqual(var1.pos_from, var2[1])
        self.assertEqual(var1.ids, var2[2])
        self.assertEqual(var1.ref, var2[3])
        self.assertEqual(var1.alt, var2[4])

    def test_read_VCF_line(self):
        with open(os.path.join(self.data_dir, "vcf_example.vcf"), "r") as vcf_file:
            vcf_handler = VCFReader(vcf_file)
            vcf_handler.read_header()
            self.assertEqual(len(vcf_handler.header.file_metadata), 7)
            self.assertEqual(len(vcf_handler.header.samples), 2)

            records = list(vcf_handler.read_records())
            self.assertEqual(len(records), 2)

            # test first record fully
            self.variant_is_equal(records[0], ("20", 9, set(), "CT", "C"))  # zero=based representation
            self.assertEqual(records[0].filters, set())
            self.assertEqual(records[0].passes_filter, True)

            self.assertEqual(len(records[0].info), 12)
            self.assertEqual(records[0].info["PP"], [3000])
            self.assertEqual(records[0].info["DP"], [250])
            self.assertEqual(records[0].info["DPR"], [140])
            self.assertEqual(records[0].info["DPF"], [110])
            self.assertEqual(records[0].info["VC"], [100])
            self.assertEqual(records[0].info["VCR"], [49])
            self.assertEqual(records[0].info["VCF"], [51])
            self.assertEqual(records[0].info["ABPV"], [0.2])
            self.assertEqual(records[0].info["SBPV"], [0.3])
            self.assertEqual(records[0].info["MQ"], [70])
            self.assertEqual(records[0].info["BR"], [31])
            self.assertEqual(records[0].info["QD"], [None])

            self.assertEqual(records[0].samples, ['sample1', 'sample2'])
            self.assertEqual(records[0].sample_info.get_field('sample1', "GT"), GenotypeCall("0/1"))
            self.assertEqual(records[0].sample_info.get_field('sample2', "GT"), GenotypeCall("1/1"))

            self.assertEqual(records[0].sample_info.get_field('sample1', 'PL'), [3000, 0, 3000])
            self.assertEqual(records[0].sample_info.get_field('sample2', 'PL'), [114, 0, 0])

            self.assertEqual(records[0].sample_info.get_field('sample1', 'GQ'), [1000])
            self.assertEqual(records[0].sample_info.get_field('sample2', 'GQ'), [None])

            # check that ordering in the dictionaries is preserved
            expected_keys = ["PP", "DP", "DPR", "DPF", "VC", "VCR",
                             "VCF", "ABPV", "SBPV", "MQ", "BR", "QD"]

            self.assertEqual(list(records[0].info.keys()), expected_keys)

            # ensure last record is still being read correctly
            self.variant_is_equal(records[-1], ("20", 10, set(), "T", "G"))

    def test_reads_simple_file(self):
        filename = os.path.join(self.work_dir, "test.vcf")

        with VCFWriterContextManager(filename) as left_vcf:
            left_vcf.write_variant(Variant("1", 1, "A", "T"))
            left_vcf.write_variant(Variant("2", 1, "A", "T"))
            left_vcf.write_variant(Variant("10", 1, "A", "T"))

        expected_variants = [
            Variant("1", 1, "A", "T"),
            Variant("2", 1, "A", "T"),
            Variant("10", 1, "A", "T"),
        ]

        with VCFReaderContextManager(filename) as vcf_reader:
            actual_variants = [record.variant for record in vcf_reader.read_records()]

        self.assertEqual(expected_variants, actual_variants)


class TestVCFStringParsing(unittest.TestCase):

    def test_should_decode_empty_VCF_string(self):
        self.assertEqual('', decode_VCF_string('""'))

    def test_should_decode_simple_VCF_string(self):
        self.assertEqual('foo', decode_VCF_string('"foo"'))

    def test_should_decode_VCF_string_with_single_double_quote(self):
        self.assertEqual('"', decode_VCF_string('"\\""'))

    def test_should_decode_VCF_string_with_single_backslash(self):
        self.assertEqual('\\', decode_VCF_string('"\\\\"'))

    def test_should_decode_complex_VCF_string(self):
        self.assertEqual(
            'abc\\def"ghi',
            decode_VCF_string('"abc\\\\def\\\"ghi"'))

    def test_should_fail_to_decode_unquoted_string(self):
        with self.assertRaisesRegex(Exception, 'expected a VCF encoded string: \'foo\''):
            print(decode_VCF_string('foo'))

    def test_should_fail_to_decode_string_with_stray_backslash(self):
        with self.assertRaisesRegex(Exception, re.escape('expected a VCF encoded string: \'"\\\\"\'')):
            print(decode_VCF_string('"\\"'))

    def test_should_fail_to_decode_string_with_unencoded_double_quote(self):
        with self.assertRaisesRegex(Exception, 'expected a VCF encoded string: \'"\""\''):
            print(decode_VCF_string('"\""'))


class TestCommaSeparatedPairParser(unittest.TestCase):

    def test_should_parse_simple_comma_separated_pairs(self):
        parsed = parse_VCF_comma_separated_pair_value('<first=foo,second=bar>')
        expected = {'first': 'foo', 'second': 'bar'}
        self.assertEqual(expected, parsed)

    def test_should_parse_empty_simple_value(self):
        parsed = parse_VCF_comma_separated_pair_value('<first=,second=bar>')
        expected = {'first': '', 'second': 'bar'}
        self.assertEqual(expected, parsed)

    def test_should_fail_to_parse_non_bracketed_string(self):
        with self.assertRaisesRegex(Exception, 'expected braced key-value pairs: \'first=foo\''):
            print(parse_VCF_comma_separated_pair_value('first=foo'))

    def test_should_parse_quoted_comma_separated_pairs(self):
        parsed = parse_VCF_comma_separated_pair_value(
            '<first="foo",second="bar">')
        expected = {'first': '"foo"', 'second': '"bar"'}
        self.assertEqual(expected, parsed)

    def test_should_parse_empty_quoted_value(self):
        parsed = parse_VCF_comma_separated_pair_value('<first="">')
        expected = {'first': '""'}
        self.assertEqual(expected, parsed)

    def test_should_parse_values_with_quoted_commas(self):
        parsed = parse_VCF_comma_separated_pair_value('<first="foo,bar">')
        expected = {'first': '"foo,bar"'}
        self.assertEqual(expected, parsed)

    def test_should_parse_values_with_quoted_double_quote(self):
        parsed = parse_VCF_comma_separated_pair_value('<first="foo\\\"bar">')
        expected = {'first': '"foo\\\"bar"'}
        self.assertEqual(expected, parsed)

    def test_should_fail_with_badly_quoted_double_quote(self):
        with self.assertRaisesRegex(Exception, 'failed to parse key-value pairs from \'<first="foo\"bar">\''):
            print(parse_VCF_comma_separated_pair_value('<first="foo\"bar">'))


class TestHeaderParsing(unittest.TestCase):

    # version parsing

    def test_should_parse_well_formatted_version(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        self.assertEqual(expected, header)

    def test_should_store_header_as_attribute_of_parser(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        self.assertEqual(header, reader.header)

    def test_should_fail_with_unexpected_version(self):
        lines = [
            '##fileformat=VCFv0.0\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(Exception, 'unexpected version: \'0.0\''):
            print(reader.read_header())

    def test_should_fail_to_parse_malformed_header_line(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '##malformed line!\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(Exception, 'failed to parse header line: \'##malformed line!\''):
            print(reader.read_header())

    def test_should_fail_if_version_is_not_defined(self):
        lines = [
            '##notFileformat=foo\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(Exception, 'unrecognised file format line: \'##notFileformat=foo\''):
            print(reader.read_header())

    # file metadata parsing

    def test_should_parse_well_formatted_file_metadata(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '##fileDate=2013-07-08\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.file_metadata['fileDate'] = '2013-07-08'
        self.assertEqual(expected, header)

    # info data parsing

    def test_should_parse_minimal_info_header_fields(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=key,Number=1,Type=String,Description="description">\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.set_info_data('key', '1', 'String', 'description')
        self.assertEqual(expected, header)

    def test_should_parse_all_info_header_fields(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '##INFO=<ID=key,Number=1,Type=String,Description="description",Source="foo",Version="bar">\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.set_info_data(
            'key',
            '1',
            'String',
            'description',
            'foo',
            'bar')
        self.assertEqual(expected, header)

    # sample data parsing

    def test_should_parse_valid_sample_header_fields(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '##FORMAT=<ID=key,Number=1,Type=String,Description="description">\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.set_sample_data('key', '1', 'String', 'description')
        self.assertEqual(expected, header)

    # filter parsing

    def test_should_parse_valid_filter_header_fields(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '##FILTER=<ID=key,Description="description">\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.set_filter('key', 'description')
        self.assertEqual(expected, header)

    # contig parsing

    def test_should_parse_valid_contig_header_fields(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '##contig=<ID=key,length=666>\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.set_contig('key', 666)
        self.assertEqual(expected, header)

    # column headers + sample names

    def test_should_parse_required_column_headers(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        self.assertEqual(expected, header)

    def test_should_fail_without_required_column_headers(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(
                Exception,
                re.escape("expected column header line: '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER'")
        ):
            print(reader.read_header())

    def test_should_parse_column_headers_with_format_but_no_samples(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        self.assertEqual(expected, header)

    def test_should_parse_column_headers_with_complex_sample_names(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOWEN_TOBY-RHYS.JONES\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.samples = ['OWEN_TOBY-RHYS.JONES']
        self.assertEqual(expected, header)

    def test_should_not_parse_column_headers_with_sample_names_containing_white_space(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOWEN JONES\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(
                Exception,
                re.escape(
                    'expected column header line: '
                    '\'#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tOWEN JONES\''
                )
        ):
            print(reader.read_header())

    def test_should_fail_with_malformed_format_column_header(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFOO\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(
                Exception,
                re.escape('expected column header line: \'#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFOO\'')
        ):
            print(reader.read_header())

    def test_should_parse_column_headers_with_samples(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFOO\tBAR\n',
        ]
        reader = VCFReader(iter(lines))

        header = reader.read_header()

        expected = Schema()
        expected.samples.append('FOO')
        expected.samples.append('BAR')
        self.assertEqual(expected, header)

    def test_should_fail_if_column_header_line_is_missing(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            'the line after the header\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(Exception, 'expected column header line: \'the line after the header\''):
            print(reader.read_header())

    def test_should_fail_on_unexpected_EOF(self):
        lines = [
            '##fileformat=VCFv4.2\n',
        ]
        reader = VCFReader(iter(lines))

        with self.assertRaisesRegex(Exception, 'unexpected EOF'):
            print(reader.read_header())


class TestRecordParsing(unittest.TestCase):

    # version parsing

    def test_should_parse_single_record(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
            'chr0\t0\t.\tP\tQ\t0\tPASS\t\n',
        ]
        reader = VCFReader(iter(lines))

        record_count = len(list(reader.read_records()))

        self.assertEqual(1, record_count)

    def test_should_parse_header_when_parsing_records(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
            'chr0\t0\t.\tP\tQ\t0\tPASS\t\n',
        ]
        reader = VCFReader(iter(lines))

        self.assertIsNone(reader.header)
        list(reader.read_records())
        self.assertIsNotNone(reader.header)

    def test_should_parse_empty_file(self):
        lines = [
            '##fileformat=VCFv4.2\n',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n',
        ]
        reader = VCFReader(iter(lines))

        record_count = len(list(reader.read_records()))

        self.assertEqual(0, record_count)
