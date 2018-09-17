# All content Copyright (C) 2018 Genomics plc
from wecall.vcfutils import schema, fieldmetadata
import unittest


# Note: using iterator-based access to decouple from object type
class TestSchemaDataModel(unittest.TestCase):

    # simple data:

    def test_should_contain_mutable_file_metadata(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(schema_data.file_metadata.items()))

        schema_data.file_metadata['key'] = 'value'
        self.assertEqual([('key', 'value')], list(
            schema_data.file_metadata.items()))
        self.assertEqual('value', schema_data.file_metadata['key'])

        del schema_data.file_metadata['key']
        self.assertEqual([], list(schema_data.file_metadata.items()))

    def test_should_contain_mutable_samples_sequence(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(iter(schema_data.samples)))

        schema_data.samples.append('sample_name')
        self.assertEqual(['sample_name'], list(iter(schema_data.samples)))

        schema_data.samples.remove('sample_name')
        self.assertEqual([], list(iter(schema_data.samples)))

    # info data

    def test_should_contain_mutable_info_data_with_required_fields(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(schema_data.iter_info_data()))

        schema_data.set_info_data(
            'key',
            number=1,
            data_type='String',
            description='description'
        )

        expected_data = fieldmetadata.InfoMetadata(
            number=1,
            data_type='String',
            description='description'
        )
        self.assertEqual([('key', expected_data)],
                         list(schema_data.iter_info_data()))
        self.assertEqual(expected_data, schema_data.get_info_data('key'))

        schema_data.del_info_data('key')
        self.assertEqual([], list(schema_data.iter_info_data()))

    def test_should_contain_mutable_info_data_with_all_fields(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(schema_data.iter_info_data()))

        schema_data.set_info_data(
            'key',
            number=1,
            data_type='String',
            description='description',
            source='source',
            version='version'
        )

        expected_data = fieldmetadata.InfoMetadata(
            number=1,
            data_type='String',
            description='description',
            source='source',
            version='version'
        )
        self.assertEqual([('key', expected_data)],
                         list(schema_data.iter_info_data()))
        self.assertEqual(expected_data, schema_data.get_info_data('key'))

        schema_data.del_info_data('key')
        self.assertEqual([], list(schema_data.iter_info_data()))

    # sample data

    def test_should_contain_mutable_sample_data(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(schema_data.iter_sample_data()))

        schema_data.set_sample_data(
            'key',
            number=1,
            data_type='String',
            description='description'
        )

        expected_data = fieldmetadata.SampleMetadata(
            number=1,
            data_type='String',
            description='description'
        )
        self.assertEqual([('key', expected_data)],
                         list(schema_data.iter_sample_data()))
        self.assertEqual(expected_data, schema_data.get_sample_data('key'))

        schema_data.del_sample_data('key')
        self.assertEqual([], list(schema_data.iter_sample_data()))

    # filters

    def test_should_contain_mutable_filter_data(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(schema_data.iter_filters()))

        schema_data.set_filter('key', description='description')

        expected_data = fieldmetadata.FilterMetadata(description='description')
        self.assertEqual([('key', expected_data)],
                         list(schema_data.iter_filters()))
        self.assertEqual(expected_data, schema_data.get_filter('key'))

        schema_data.del_filter('key')
        self.assertEqual([], list(schema_data.iter_filters()))

    # contigs

    def test_should_contain_mutable_contig_data_with_required_fields(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(schema_data.iter_contigs()))

        schema_data.set_contig('key')

        expected_data = fieldmetadata.ContigMetadata()
        self.assertEqual([('key', expected_data)],
                         list(schema_data.iter_contigs()))
        self.assertEqual(expected_data, schema_data.get_contig('key'))

        schema_data.del_contig('key')
        self.assertEqual([], list(schema_data.iter_contigs()))

    def test_should_contain_mutable_contig_data_with_all_fields(self):
        schema_data = schema.Schema()
        self.assertEqual([], list(schema_data.iter_contigs()))

        schema_data.set_contig('key', length=100000)

        expected_data = fieldmetadata.ContigMetadata(length=100000)
        self.assertEqual([('key', expected_data)],
                         list(schema_data.iter_contigs()))
        self.assertEqual(expected_data, schema_data.get_contig('key'))

        schema_data.del_contig('key')
        self.assertEqual([], list(schema_data.iter_contigs()))
