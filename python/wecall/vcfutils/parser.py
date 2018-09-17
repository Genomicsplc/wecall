# All content Copyright (C) 2018 Genomics plc

import datetime
import gzip
import re
from wecall.genomics.chromosome import chromosome_comp
from wecall.utils.tabix_wrapper import TabixWrapper
import os
import wecall.vcfutils.schema
import wecall.vcfutils.record
from wecall.utils.interval import WHOLE_GENOME


ADAPTER_KEY = 'ADAPTER'


def read_records_functor(vcf_handler):
    return vcf_handler.read_records


def read_variants_functor(vcf_handler):

    def read_variants(*args, **kwargs):
        return (
            record.variant for record in vcf_handler.read_records(
                *args, **kwargs))

    return read_variants


class VCFReader(object):

    def __init__(self, input_stream):
        self.header = None
        self.__input_stream = input_stream
        self.__next_line = None  # one line of lookahead
        self.__lines = self.__readline()

    def read_header(self):
        header = wecall.vcfutils.schema.Schema()

        # get the version from a VCF file
        line = next(self.__lines).strip()
        version_match = re.match(
            '^##fileformat=VCFv(?P<version>[\d\.]+)$', line)
        if version_match:
            version = version_match.group('version')
        else:
            raise Exception('unrecognised file format line: {!r}'.format(line))

        if version not in {
            '4.0',
            '4.1',
                '4.2'}:  # replace this test with a layer of indirection
            raise Exception('unexpected version: {!r}'.format(version))

        header.set_vcf_format(version)

        # parse a VCF4.2 format header:
        header_regex = re.compile('^##(?P<key>[-_~\w\d.]+)=(?P<value>.*)$')
        while True:
            line = next(self.__lines, None)
            if line is None:
                break
            line = line.strip()
            header_match = header_regex.match(line)
            if header_match:
                key, value = header_match.group('key', 'value')
                if key == 'INFO':  # replace if-elif-else with layer of indirection
                    data = parse_VCF_comma_separated_pair_value(value)
                    source = decode_VCF_string(
                        data['Source']) if 'Source' in data else None
                    version = decode_VCF_string(
                        data['Version']) if 'Version' in data else None
                    header.set_info_data(
                        data['ID'],
                        data['Number'],
                        data['Type'],
                        decode_VCF_string(data['Description']),
                        source,
                        version,
                    )
                elif key == 'FORMAT':
                    data = parse_VCF_comma_separated_pair_value(value)
                    header.set_sample_data(
                        data['ID'],
                        data['Number'],
                        data['Type'],
                        decode_VCF_string(data['Description']),
                    )
                elif key == 'FILTER':
                    data = parse_VCF_comma_separated_pair_value(value)
                    header.set_filter(
                        data['ID'],
                        decode_VCF_string(data['Description']),
                    )
                elif key == 'contig':
                    data = parse_VCF_comma_separated_pair_value(value)
                    header.set_contig(
                        data['ID'],
                        get_int_or_none(data, 'length'),
                    )
                elif key == ADAPTER_KEY:
                    data = parse_VCF_comma_separated_pair_value(value)
                    if 'ID' in data:
                        try:
                            date = datetime.datetime.strptime(
                                data['date'], '%Y-%m-%dT%H:%M:%S')
                        except ValueError:
                            date = datetime.datetime.strptime(
                                data['date'], '%Y-%m-%d')
                        header.set_adapter(data['ID'], data['hash'], date)
                    else:
                        header.set_adapter(
                            data['adapters'],
                            data['githash'],
                            datetime.datetime.strptime(
                                data['date'],
                                '%Y-%m-%d'),
                        )
                else:  # general file metadata
                    header.file_metadata[key] = value
            elif line.startswith('##'):
                raise Exception(
                    'failed to parse header line: {!r}'.format(line))
            else:
                break

        if line is None:
            raise Exception('unexpected EOF')

        column_match = re.match(
            '^#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO(?:\tFORMAT(?P<samples>(?:\t\S+)*))?$',
            line)
        if column_match:
            sample_group = column_match.group('samples')
            if sample_group is not None:
                samples = [_f for _f in sample_group.split('\t') if _f]
                header.samples.extend(samples)
        else:
            raise Exception('expected column header line: {!r}'.format(line))

        self.header = header
        return header

    def read_records(self):
        if self.header is None:
            self.header = self.read_header()
        for line_number, line in enumerate(self.__readline()):
            try:
                if line is None:
                    raise StopIteration()
                for rec in wecall.vcfutils.record.read_records(
                        self.header, line):
                    yield rec
            except Exception as e:
                raise type(e)(
                    "Error processing VCF data line {}: {}".format(
                        line_number + 1, e))

    def __readline(self):
        for line in self.__input_stream:
            prevline = self.__next_line
            self.__next_line = line
            if prevline:
                yield prevline
        line, self.__next_line = self.__next_line, None
        yield line


def get_int_or_none(map, key):
    try:
        return int(map[key])
    except KeyError:
        return None


def parse_VCF_comma_separated_pair_value(value):
    braced_value_regex = re.compile('^<(?P<value>.*)>$')
    # value ::= quoted value | non-quoted value
    # quoted value ::= '"' ( escaped character | not('"') ) * '"'
    # non-quoted value ::= ( not('"', ',') not(',') * ) ?
    comma_separated_pair_regex = re.compile(
        '^(?P<key>[^=]+)=(?P<value>(?:"(?:(?:\\\\.)|[^"])*")|(?:[^,"][^,]*)?),?')
    info_match = braced_value_regex.match(value)
    if info_match:
        stripped_value = info_match.group('value')
        data = {}
        end = 0
        while end != len(stripped_value):
            csp_match = comma_separated_pair_regex.match(stripped_value[end:])
            if csp_match:
                end += len(csp_match.group(0))
                data[csp_match.group('key')] = csp_match.group('value')
            else:
                raise ValueError(
                    'failed to parse key-value pairs from {!r}'.format(value))
        return data
    else:
        raise ValueError('expected braced key-value pairs: {!r}'.format(value))


def decode_VCF_string(string):
    if not re.match('^"(?:[^"\\\\]|\\\\"|\\\\\\\\)*"$', string):
        raise ValueError('expected a VCF encoded string: {!r}'.format(string))
    else:
        return string[1:-1].replace('\\"', '"').replace('\\\\', '\\')


class RecordAndKey(object):

    def __init__(self, record, key):
        self.key = key
        self.record = record

    def __repr__(self):
        return "<Key={}, Record={!r}".format(self.key, self.record.variant)

    def __eq__(self, other):
        return self.record == other.record and self.key == other.key

    def __lt__(self, other):
        if self.record == other.record:
            return self.key < other.key
        else:
            return self.record < other.record

    def __gt__(self, other):
        if self.record == other.record:
            return self.key > other.key
        else:
            return self.record > other.record


class VCFReaderContextManager(object):

    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        if os.path.splitext(self.filename)[1] == ".gz":
            self.fp = gzip.open(self.filename, "rt")
        else:
            self.fp = open(self.filename, "r")
        self.vcf_reader = VCFReader(self.fp)
        return self.vcf_reader

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fp.close()
