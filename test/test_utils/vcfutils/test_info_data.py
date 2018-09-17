# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase

from wecall.vcfutils.record import generate_records
from wecall.vcfutils.schema import Schema

import testfixtures


class TestMalformedDeferredInfoDataParsing(TestCase):

    @testfixtures.log_capture()
    def test_should_warn_about_unrecognised_key_in_monoallelic_line(self, log):
        records = list(generate_records(Schema(), [
            'chrZ', '200', '.', 'C', 'T', '.', 'PASS', 'NEW_KEY=value'
        ]))
        for index, record in enumerate(records):
            self.assertEqual(
                (index, ['value']), (index, record.info['NEW_KEY']))
        log.check(
            ('root',
             'WARNING',
             'info field {!r} not defined in schema'.format('NEW_KEY')),
        )

    def test_should_add_default_parsing_rule_for_unknown_key_in_monoallelic_line(self):
        schema = Schema()
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'T', '.', 'PASS', 'NEW_KEY=value'
        ]))

        self.assertEqual(0, len(list(schema.iter_info_data())))
        for index, record in enumerate(records):
            self.assertEqual(
                (index, ['value']), (index, record.info['NEW_KEY']))
        self.assertEqual(1, len(list(schema.iter_info_data())))

        info_metadata = schema.get_info_data('NEW_KEY')
        self.assertEqual('.', info_metadata.number)
        self.assertEqual('String', info_metadata.data_type)
        self.assertEqual(
            'Inferred from file content during parsing',
            info_metadata.description)
        self.assertEqual('vcfutils', info_metadata.source)
        self.assertEqual('undefined', info_metadata.version)

    @testfixtures.log_capture()
    def test_should_warn_about_unrecognised_key_in_multiallelic_line(
            self,
            log):
        records = list(generate_records(Schema(), [
            'chrZ', '200', '.', 'C', 'A,T', '.', 'PASS', 'NEW_KEY=value'
        ]))
        for index, record in enumerate(records):
            self.assertEqual(
                (index, ['value']), (index, record.info['NEW_KEY']))
        log.check(
            ('root',
             'WARNING',
             'info field {!r} not defined in schema'.format('NEW_KEY')),
        )

    def test_should_add_default_parsing_rule_for_unknown_key_in_multiallelic_line(self):
        schema = Schema()
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A,T', '.', 'PASS', 'NEW_KEY=value'
        ]))

        self.assertEqual(0, len(list(schema.iter_info_data())))
        for index, record in enumerate(records):
            self.assertEqual(
                (index, ['value']), (index, record.info['NEW_KEY']))
        self.assertEqual(1, len(list(schema.iter_info_data())))

        info_metadata = schema.get_info_data('NEW_KEY')
        self.assertEqual('.', info_metadata.number)
        self.assertEqual('String', info_metadata.data_type)
        self.assertEqual(
            'Inferred from file content during parsing',
            info_metadata.description)
        self.assertEqual('vcfutils', info_metadata.source)
        self.assertEqual('undefined', info_metadata.version)

    @testfixtures.log_capture()
    def test_should_warn_about_too_few_alts_in_field_of_allelic_cardinality(
            self,
            log):
        schema = Schema()
        schema.set_info_data('key', 'A', 'String', '')
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A,T', '.', 'PASS', 'key=a'
        ]))
        expected = [['a'], [None]]
        for index, record in enumerate(records):
            self.assertEqual(expected[index], record.info['key'])
        log.check(('wecall.vcfutils.fieldmetadata', 'WARNING',
                   'expected 2 items in {!r}'.format([['a']])), )

    @testfixtures.log_capture()
    def test_should_warn_about_too_many_alts_in_field_of_allelic_cardinality(
            self,
            log):
        schema = Schema()
        schema.set_info_data('key', 'A', 'String', '')
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A,T', '.', 'PASS', 'key=a,b,c'
        ]))
        expected = [['a'], ['b']]
        for index, record in enumerate(records):
            self.assertEqual(expected[index], record.info['key'])
        log.check(('wecall.vcfutils.fieldmetadata', 'WARNING',
                   'expected 2 items in {!r}'.format([['a'], ['b'], ['c']])), )
