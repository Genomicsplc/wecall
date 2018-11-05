# All content Copyright (C) 2018 Genomics plc
import unittest

from wecall.common.exceptions import weCallException
from wecall.vcfutils.record import generate_records
from wecall.vcfutils.sample_data import SampleData
from wecall.vcfutils.genotype_call import GenotypeCall
from wecall.vcfutils.schema import Schema

import testfixtures


class TestSampleData(unittest.TestCase):

    def test_default_field_value_is_assigned_when_sample_data_is_constructed(self):
        sample_data = SampleData(['key1'], ['sample_name'])
        self.assertEqual(sample_data.get_field('sample_name', 'key1'), [])

    def test_genotype_field_default_value_is_assigned_when_sample_data_is_constructed(self):
        sample_data = SampleData(['GT'], ['sample_name'])
        self.assertEqual(
            sample_data.get_field(
                'sample_name',
                'GT'),
            GenotypeCall("./."))

    def test_default_values_are_assigned_when_sample_data_is_constructed(self):
        sample_data = SampleData(['GT', 'key1', 'key2'], [
                                 'sample_name1', 'sample_name2'])
        self.assertEqual(
            sample_data.get_field(
                'sample_name1',
                'GT'),
            GenotypeCall("./."))
        self.assertEqual(
            sample_data.get_field(
                'sample_name2',
                'GT'),
            GenotypeCall("./."))
        self.assertEqual(sample_data.get_field('sample_name1', 'key1'), [])
        self.assertEqual(sample_data.get_field('sample_name2', 'key1'), [])
        self.assertEqual(sample_data.get_field('sample_name1', 'key2'), [])
        self.assertEqual(sample_data.get_field('sample_name2', 'key2'), [])

    def test_has_sample_reports_expected_value(self):
        sample_data = SampleData(['key1'], ['sample_name'])
        self.assertTrue(sample_data.has_sample('sample_name'))
        self.assertFalse(sample_data.has_sample('missing_sample_name'))

    def test_has_genotype_key_should_report_expected_value(self):
        sample_data = SampleData(['genotype_key'], ['sample_name'])
        self.assertTrue(sample_data.has_genotype_key('genotype_key'))
        self.assertFalse(sample_data.has_genotype_key('missing_genotype_key'))

    def test_has_genotype_keys_should_support_multiple_keys(self):
        sample_data = SampleData(
            ['genotype_key1', 'genotype_key2'], ['sample_name'])
        self.assertTrue(sample_data.has_genotype_key('genotype_key1'))
        self.assertTrue(sample_data.has_genotype_key('genotype_key2'))

    def test_should_add_sample_data(self):
        sample_data = SampleData(['genotype_key1'], ['sample_name'])
        sample_data.add_sample_data('sample_name', 'genotype_key1', [1])
        self.assertEqual(
            sample_data.get_field(
                'sample_name',
                'genotype_key1'),
            [1])

    def test_should_allow_multiple_calls_to_add_sample_data(self):
        sample_data = SampleData(
            ['genotype_key1', 'genotype_key2'], ['sample_name'])
        sample_data.add_sample_data('sample_name', 'genotype_key1', [1])
        sample_data.add_sample_data('sample_name', 'genotype_key2', [2])
        sample_data.add_sample_data('sample_name', 'genotype_key1', [3, 4])
        self.assertEqual(
            sample_data.get_field(
                'sample_name', 'genotype_key1'), [
                3, 4])
        self.assertEqual(
            sample_data.get_field(
                'sample_name',
                'genotype_key2'),
            [2])

    def test_should_allow_multiple_samples_for_add_sample_data(self):
        sample_data = SampleData(
            ['genotype_key1'], [
                'sample_name1', 'sample_name2'])
        sample_data.add_sample_data('sample_name1', 'genotype_key1', [1])
        sample_data.add_sample_data('sample_name2', 'genotype_key1', [3, 4])
        self.assertEqual(
            sample_data.get_field(
                'sample_name1',
                'genotype_key1'),
            [1])
        self.assertEqual(
            sample_data.get_field(
                'sample_name2', 'genotype_key1'), [
                3, 4])

    def test_should_raise_when_adding_sample_data_to_missing_key(self):
        sample_data = SampleData(['key'], ['sample_name'])
        self.assertRaisesRegex(
            weCallException,
            "Missing key missing_key when adding sample data.",
            sample_data.add_sample_data,
            'sample_name',
            'missing_key',
            [1],
        )

    def test_should_raise_when_adding_sample_data_to_missing_sample(self):
        sample_data = SampleData(['key'], ['sample_name'])
        self.assertRaisesRegex(
            weCallException,
            "Missing sample name missing_sample_name supplied when adding sample data.",
            sample_data.add_sample_data,
            'missing_sample_name',
            'key',
            [1],
        )

    def test_should_raise_when_adding_wrong_genotype_data(self):
        sample_data = SampleData(['GT'], ['sample_name'])
        self.assertRaisesRegex(
            weCallException,
            "Genotype field must be a GenotypeCall.",
            sample_data.add_sample_data,
            'sample_name',
            'GT',
            [1],
        )


class TestSampleDataGetGenotypeLikelihood(unittest.TestCase):
    def test_gets_exact_values_if_key_is_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.add_sample_data(sample_name, 'GL', [-0.1, -0.2, -0.3])
        self.assertEqual(sample_data.get_genotype_likelihoods(
            sample_name), [-0.1, -0.2, -0.3])

    def test_gets_exact_values_if_key_is_PL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['PL'], [sample_name])
        sample_data.add_sample_data(sample_name, 'PL', [1, 2, 3])
        self.assertEqual(sample_data.get_genotype_likelihoods(
            sample_name), [-0.1, -0.2, -0.3])

    def test_gets_dot_if_key_is_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.add_sample_data(sample_name, 'GL', '.')
        self.assertEqual(
            sample_data.get_genotype_likelihoods(sample_name), '.')

    def test_gets_list_of_dot_if_key_is_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.add_sample_data(sample_name, 'GL', ['.', '.', '.'])
        self.assertEqual(
            sample_data.get_genotype_likelihoods(sample_name), [
                None, None, None])

    def test_gets_none_if_key_is_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.add_sample_data(sample_name, 'GL', None)
        self.assertEqual(
            sample_data.get_genotype_likelihoods(sample_name), None)

    def test_gets_list_of_none_if_key_is_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.add_sample_data(sample_name, 'GL', [None, None, None])
        self.assertEqual(
            sample_data.get_genotype_likelihoods(sample_name), [
                None, None, None])

    def test_gets_dot_if_key_is_PL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['PL'], [sample_name])
        sample_data.add_sample_data(sample_name, 'PL', '.')
        self.assertEqual(
            sample_data.get_genotype_likelihoods(sample_name), '.')

    def test_sample_data_copes_with_mixed_dot_missing_values_in_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.add_sample_data(
            sample_name, 'GL', [-0.1, '.', -0.2, None, -0.3])
        self.assertEqual(sample_data.get_genotype_likelihoods(
            sample_name), [-0.1, None, -0.2, None, -0.3])

    def test_sample_data_copes_with_mixed_missing_values_in_PL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['PL'], [sample_name])
        sample_data.add_sample_data(
            sample_name, 'PL', [-0.1, '.', -0.2, None, -0.3])
        self.assertEqual(sample_data.get_genotype_likelihoods(
            sample_name), [0.01, None, 0.02, None, 0.03])


class TestSampleDataSetGenotypeLikelihood(unittest.TestCase):
    def test_gets_exact_values_if_key_is_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.set_genotype_likelihoods(sample_name, [-0.1, -0.2, -0.3])
        self.assertEqual(sample_data.get_field(
            sample_name, 'GL'), [-0.1, -0.2, -0.3])

    def test_gets_exact_values_if_key_is_PL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['PL'], [sample_name])
        sample_data.set_genotype_likelihoods(sample_name, [-0.1, -0.2, -0.3])
        self.assertEqual(sample_data.get_field(sample_name, 'PL'), [1, 2, 3])

    def test_gets_dot_if_key_is_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.set_genotype_likelihoods(sample_name, '.')
        self.assertEqual(sample_data.get_field(sample_name, 'GL'), '.')

    def test_gets_dot_if_key_is_PL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['PL'], [sample_name])
        sample_data.set_genotype_likelihoods(sample_name, '.')
        self.assertEqual(sample_data.get_field(sample_name, 'PL'), '.')

    def test_sample_data_copes_with_mixed_dot_missing_values_in_GL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GL'], [sample_name])
        sample_data.set_genotype_likelihoods(
            sample_name, [-0.1, '.', -0.2, None, -0.3])
        self.assertEqual(sample_data.get_field(
            sample_name, 'GL'), [-0.1, None, -0.2, None, -0.3])

    def test_sample_data_copes_with_mixed_missing_values_in_PL(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['PL'], [sample_name])
        sample_data.set_genotype_likelihoods(
            sample_name, [-0.1, '.', -0.2, None, -0.3])
        self.assertEqual(
            sample_data.get_field(
                sample_name, 'PL'), [
                1.0, None, 2.0, None, 3.0])


class TestSampleDataGetReadDepth(unittest.TestCase):
    def test_gets_exact_values_if_key_is_DP(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['DP'], [sample_name])
        sample_data.add_sample_data(sample_name, 'DP', [100])
        self.assertEqual(sample_data.get_read_depth(sample_name), [100])

    def test_gets_exact_values_if_key_is_NR(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['NR'], [sample_name])
        sample_data.add_sample_data(sample_name, 'NR', [100])
        self.assertEqual(sample_data.get_read_depth(sample_name), [100])


class TestSampleDataGetVariantSupport(unittest.TestCase):
    def test_gets_exact_values_if_key_is_AD(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['AD'], [sample_name])
        sample_data.add_sample_data(sample_name, 'AD', [None, 100])
        self.assertEqual(sample_data.get_variant_support(sample_name), [100])

    def test_gets_exact_values_if_key_is_NV(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['NV'], [sample_name])
        sample_data.add_sample_data(sample_name, 'NV', [100])
        self.assertEqual(sample_data.get_variant_support(sample_name), [100])


class TestSampleDataGetGenotypeQuality(unittest.TestCase):
    def test_gets_value_for_GQ_key(self):
        sample_name = 'sample_name'
        sample_data = SampleData(['GQ'], [sample_name])
        sample_data.add_sample_data(sample_name, 'GQ', [2.3])
        self.assertEqual(sample_data.get_genotype_quality(sample_name), [2.3])


class TestGenotypeDataView(unittest.TestCase):

    def setUp(self):
        self.sample_data = SampleData(
            ['GT', 'key'], ['sample_name1', 'sample_name2'])
        self.sample_data.add_sample_data("sample_name1", "key", [1, 2])
        self.sample_data.add_sample_data(
            "sample_name2", "GT", GenotypeCall("0/1"))

    def test_contains_method_returns_expected_value_sample1(self):
        genotype_data = self.sample_data.get_genotype_data("sample_name1")
        self.assertNotIn("cheesecake", genotype_data)
        self.assertNotIn("sample_name1", genotype_data)
        self.assertIn("GT", genotype_data)
        self.assertIn("key", genotype_data)

    def test_contains_method_returns_expected_value_sample2(self):
        genotype_data = self.sample_data.get_genotype_data("sample_name2")
        self.assertIn("GT", genotype_data)
        self.assertIn("key", genotype_data)

    def test_getitem_method_returns_expected_value(self):
        genotype_data = self.sample_data.get_genotype_data("sample_name1")
        self.assertEqual(genotype_data["GT"], GenotypeCall("./."))
        self.assertEqual(genotype_data["key"], [1, 2])
        genotype_data = self.sample_data.get_genotype_data("sample_name2")
        self.assertEqual(genotype_data["GT"], GenotypeCall("0/1"))
        self.assertEqual(genotype_data["key"], [])

    def test_keys_method_returns_expected_data(self):
        genotype_data = self.sample_data.get_genotype_data("sample_name1")
        self.assertEqual(list(genotype_data.keys()), ["GT", "key"])
        genotype_data = self.sample_data.get_genotype_data("sample_name2")
        self.assertEqual(list(genotype_data.keys()), ["GT", "key"])

    def test_values_method_returns_expected_data(self):
        genotype_data = self.sample_data.get_genotype_data("sample_name1")
        self.assertEqual(list(genotype_data.values()),
                         [GenotypeCall("./."), [1, 2]])
        genotype_data = self.sample_data.get_genotype_data("sample_name2")
        self.assertEqual(list(genotype_data.values()),
                         [GenotypeCall("0/1"), []])


class TestGenerateRecords(unittest.TestCase):

    def test_should_split_genotype_likelihood_properly(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('GL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A,T', '.', 'PASS', '.', 'GT:GL', '0/1:1,2,3,4,5,6'
        ]))

        self.assertEqual(
            GenotypeCall("0/1"),
            records[0].sample_info.get_field(
                'foo',
                'GT'))
        self.assertEqual(
            [1.0, 2.0, 3.0], records[0].sample_info.get_field('foo', 'GL'))
        self.assertEqual(
            GenotypeCall("0/0"),
            records[1].sample_info.get_field(
                'foo',
                'GT'))
        self.assertEqual(
            [1.0, 4.0, 6.0], records[1].sample_info.get_field('foo', 'GL'))

    def test_should_drop_genotype_likelihood_with_mismatch_ploidy(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('GL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A,T', '.', 'PASS', '.', 'GT:GL', '0/1:1,2,3,4'
        ]))

        self.assertEqual(
            GenotypeCall("0/1"),
            records[0].sample_info.get_field(
                'foo',
                'GT'))
        self.assertEqual([None, None, None],
                         records[0].sample_info.get_field('foo', 'GL'))
        self.assertEqual(
            GenotypeCall("0/0"),
            records[1].sample_info.get_field(
                'foo',
                'GT'))
        self.assertEqual([None, None, None],
                         records[1].sample_info.get_field('foo', 'GL'))

    @testfixtures.log_capture()
    def test_should_warn_when_GT_is_not_present(self, log):
        schema = Schema()
        schema.set_sample_data('GL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A,T', '.', 'PASS', '.', 'GL', '1,2,3'
        ]))
        for index, record in enumerate(records):
            self.assertEqual(
                (index, [
                    '1', '2', '3']), (index, record.sample_info.get_field(
                        'foo', 'GL')))
        log.check(('wecall.vcfutils.fieldmetadata', 'WARNING',
                   'Unknown ploidy when parsing genotype likelihood'), )


class TestHasNoLikelihoods(unittest.TestCase):
    def test_should_return_true_if_no_GL_or_PL_present(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A', '.', 'PASS', '.', 'GT', '0/1'
        ]))
        self.assertTrue(records[0].sample_info.has_no_likelihoods())

    def test_should_return_true_if_all_likelihoods_are_none_for_GL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('GL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A', '.', 'PASS', '.', 'GT:GL', '0/1:.,.,.'
        ]))
        self.assertTrue(records[0].sample_info.has_no_likelihoods())

    def test_should_return_true_if_all_likelihoods_are_none_for_PL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('PL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A', '.', 'PASS', '.', 'GT:PL', '0/1:.,.,.'
        ]))
        self.assertTrue(records[0].sample_info.has_no_likelihoods())

    def test_should_return_false_if_all_likelihoods_is_not_none_for_PL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('PL', 'G', 'Float', '')
        schema.samples = ['foo']

        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A', '.', 'PASS', '.', 'GT:PL', '0/1:90,10,12'
        ]))
        self.assertFalse(records[0].sample_info.has_no_likelihoods())

    def test_should_return_false_if_all_likelihoods_is_not_none_for_GL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('GL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A', '.', 'PASS', '.', 'GT:GL', '0/1:90,1,120'
        ]))
        self.assertFalse(records[0].sample_info.has_no_likelihoods())

    def test_should_return_false_if_any_likelihoods_is_not_none_for_PL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('PL', 'G', 'Float', '')
        schema.samples = ['foo']

        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A', '.', 'PASS', '.', 'GT:PL', '0/1:90,.,.'
        ]))
        self.assertFalse(records[0].sample_info.has_no_likelihoods())

    def test_should_return_false_if_any_likelihoods_is_not_none_for_GL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('GL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema, [
            'chrZ', '200', '.', 'C', 'A', '.', 'PASS', '.', 'GT:GL', '0/1:90,.,.'
        ]))
        self.assertFalse(records[0].sample_info.has_no_likelihoods())

    def test_should_return_false_if_one_sample_okay_for_PL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('PL', 'G', 'Float', '')
        schema.samples = ['foo', 'bar']
        records = list(generate_records(schema,
                                        ['chrZ',
                                         '200',
                                         '.',
                                         'C',
                                         'A',
                                         '.',
                                         'PASS',
                                         '.',
                                         'GT:PL',
                                         '0/1:90,1,120',
                                         '0/1:.,.,.']))
        self.assertFalse(records[0].sample_info.has_no_likelihoods())

    def test_should_return_false_if_one_sample_okay_for_GL(self):
        schema = Schema()
        schema.set_sample_data('GT', '1', 'String', '')
        schema.set_sample_data('GL', 'G', 'Float', '')
        schema.samples = ['foo']
        records = list(generate_records(schema,
                                        ['chrZ',
                                         '200',
                                         '.',
                                         'C',
                                         'A',
                                         '.',
                                         'PASS',
                                         '.',
                                         'GT:GL',
                                         '0/1:90,1,120',
                                         '0/1:.,.,.']))
        self.assertFalse(records[0].sample_info.has_no_likelihoods())
