# All content Copyright (C) 2018 Genomics plc


import datetime
from wecall.vcfutils.record import vcf_row_from_record, vcf_row_from_variant
from wecall.vcfutils.schema import Schema
from wecall.vcfutils.parser import ADAPTER_KEY


def write_record_functor(vcf_handler):
    return vcf_handler.write_record


def write_variant_functor(vcf_handler):
    return vcf_handler.write_variant


class VCFWriter(object):

    def __init__(self, stream):
        self.__fp = stream

    def write_header(self, header):
        # TODO: handle format version number in a sensible way
        print(
            '##fileformat=VCFv{}'.format(
                "4.2" if header.vcf_format is None else header.vcf_format),
            file=self.__fp)
        for key, value in list(header.file_metadata.items()):
            print(
                '##{key!s}={value!s}'.format(
                    key=key,
                    value=value),
                file=self.__fp)
        for key, value in header.iter_info_data():
            info_items = [
                'ID={}'.format(key),
                'Number={}'.format(value.number),
                'Type={}'.format(value.data_type),
                'Description={}'.format(encode_VCF_string(value.description)),
            ]
            if value.source is not None:
                info_items.append(
                    'Source={}'.format(
                        encode_VCF_string(
                            value.source)))
            if value.source is not None:
                info_items.append(
                    'Version={}'.format(
                        encode_VCF_string(
                            value.version)))
            print('##INFO=<{}>'.format(','.join(info_items)), file=self.__fp)
        for key, value in header.iter_sample_data():
            print(
                '##FORMAT=<ID={},Number={},Type={},Description={}>'.format(
                    key,
                    value.number,
                    value.data_type,
                    encode_VCF_string(
                        value.description)),
                file=self.__fp)
        for key, value in header.iter_filters():
            print('##FILTER=<ID={},Description={}>'.format(
                key, encode_VCF_string(value.description)
            ), file=self.__fp)
        for key, value in header.iter_contigs():
            print('##contig=<ID={},length={}>'.format(
                key, value.length
            ), file=self.__fp)
        for adapter in header.iter_adapters():
            print(
                '##{}=<ID={},date={},hash={}>'.format(
                    ADAPTER_KEY,
                    adapter.adapter,
                    adapter.date.strftime('%F'),
                    adapter.hash),
                file=self.__fp)
        columns = [
            'CHROM',
            'POS',
            'ID',
            'REF',
            'ALT',
            'QUAL',
            'FILTER',
            'INFO']
        if header.samples:
            columns.append('FORMAT')
            columns.extend(header.samples)
        print('#' + '\t'.join(columns), file=self.__fp)

    def write_variant(self, var, *args, **kwargs):
        self.__write(vcf_row_from_variant, var, *args, **kwargs)

    def write_record(self, rec, *args, **kwargs):
        self.__write(vcf_row_from_record, rec, *args, **kwargs)

    def write_variants(self, record_stream):
        for var in record_stream:
            self.write_variant(var)

    def write_records(self, record_stream):
        for rec in record_stream:
            self.write_record(rec)

    def __write(self, rec_type_formatter, rec_ob, *args, **kwargs):
        self.__fp.write(
            "{}\n".format(
                rec_type_formatter(
                    rec_ob,
                    *args,
                    **kwargs)))


def encode_VCF_string(string):
    return '"' + string.replace('\\', '\\\\').replace('"', '\\"') + '"'


class VCFWriterContextManager(object):

    def __init__(self, filename, header=None):
        self.filename = filename
        self.header = header

    def __enter__(self):
        self.fp = open(self.filename, 'w')
        if self.header is None:
            self.header = Schema()
            self.header.file_metadata['fileDate'] = datetime.date.today(
            ).strftime('%F')
        self.vcf_writer = VCFWriter(self.fp)
        self.vcf_writer.write_header(self.header)
        return self.vcf_writer

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fp.close()
