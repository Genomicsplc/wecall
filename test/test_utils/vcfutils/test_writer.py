# All content Copyright (C) 2018 Genomics plc
from io import StringIO
import unittest
import datetime
from wecall.vcfutils.schema import Schema
from wecall.vcfutils.writer import encode_VCF_string, VCFWriter


class TestVCFWriter(unittest.TestCase):

    def test_should_write_empty_file_containing_expected_version_number(self):
        mock_file = StringIO()
        empty_schema = Schema()
        writer = VCFWriter(mock_file)
        writer.write_header(empty_schema)
        expected_file = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        self.assertEqual(expected_file, mock_file.getvalue())

    def test_should_write_file_metadata_in_expected_format(self):
        mock_file = StringIO()
        date = datetime.datetime.utcnow().strftime('%F')
        schema = Schema()
        schema.file_metadata['fileDate'] = date

        writer = VCFWriter(mock_file)
        writer.write_header(schema)

        expected_file = '##fileformat=VCFv4.2\n' \
                        '##fileDate={date!s}\n' \
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
            .format(date=date)
        self.assertEqual(expected_file, mock_file.getvalue())

    def test_should_write_info_data_in_expected_format(self):
        mock_file = StringIO()
        schema = Schema()
        schema.set_info_data('key', '1', 'String', 'sample info field')

        writer = VCFWriter(mock_file)
        writer.write_header(schema)

        expected_file = '##fileformat=VCFv4.2\n' \
                        '##INFO=<ID=key,Number=1,Type=String,Description="sample info field">\n' \
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        self.assertEqual(expected_file, mock_file.getvalue())

    def test_should_write_sample_data_in_expected_format(self):
        mock_file = StringIO()
        schema = Schema()
        schema.set_sample_data('key', '1', 'String', 'a sample field')

        writer = VCFWriter(mock_file)
        writer.write_header(schema)

        expected_file = '##fileformat=VCFv4.2\n' \
                        '##FORMAT=<ID=key,Number=1,Type=String,Description="a sample field">\n' \
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        self.assertEqual(expected_file, mock_file.getvalue())

    def test_should_write_filter_in_expected_format(self):
        mock_file = StringIO()
        schema = Schema()
        schema.set_filter('key', 'a filter')

        writer = VCFWriter(mock_file)
        writer.write_header(schema)

        expected_file = '##fileformat=VCFv4.2\n' \
                        '##FILTER=<ID=key,Description="a filter">\n' \
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        self.assertEqual(expected_file, mock_file.getvalue())

    def test_should_write_contig_in_expected_format(self):
        mock_file = StringIO()
        schema = Schema()
        schema.set_contig('key', 666)

        writer = VCFWriter(mock_file)
        writer.write_header(schema)

        expected_file = '##fileformat=VCFv4.2\n' \
                        '##contig=<ID=key,length=666>\n' \
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        self.assertEqual(expected_file, mock_file.getvalue())

    def test_should_write_sample_names_in_column_header_line(self):
        mock_file = StringIO()
        schema = Schema()
        schema.samples.append('FOO')

        writer = VCFWriter(mock_file)
        writer.write_header(schema)

        expected_file = '##fileformat=VCFv4.2\n' \
                        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFOO\n'
        self.assertEqual(expected_file, mock_file.getvalue())


class TestVCFStringWriting(unittest.TestCase):

    def test_should_encode_empty_VCF_string(self):
        self.assertEqual('""', encode_VCF_string(''))

    def test_should_encode_simple_VCF_string(self):
        self.assertEqual('"foo"', encode_VCF_string('foo'))

    def test_should_encode_VCF_string_with_single_double_quote(self):
        self.assertEqual('"\\""', encode_VCF_string('"'))

    def test_should_encode_VCF_string_with_single_backslash(self):
        self.assertEqual('"\\\\"', encode_VCF_string('\\'))

    def test_should_encode_complex_VCF_string(self):
        self.assertEqual(
            '"abc\\\\def\\\"ghi"',
            encode_VCF_string('abc\\def"ghi'))
