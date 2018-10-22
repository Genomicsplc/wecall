# All content Copyright (C) 2018 Genomics plc
import unittest

import testfixtures

from wecall.vcfutils import fieldmetadata
from wecall.common.exceptions import weCallException
from wecall.vcfutils.fieldmetadata import make_split_sample_alt_func
from wecall.vcfutils.genotype_call import GenotypeCall


class VCFFieldMetadataTest(unittest.TestCase):

    def test_parse_float_field(self):
        field = fieldmetadata.FieldMetadata(
            'PP', 'A', 'Float', 'Posterior probability')
        self.assertEqual(field.name, "PP")
        self.assertEqual(field.description, "Posterior probability")
        self.assertEqual(field.source, None)
        self.assertEqual(field.version, None)
        self.assertEqual(field.extract_data("12.3", 1, None), [[12.3]])
        self.assertEqual(
            field.extract_data(
                "12.3,78.9", 2, None), [
                [12.3], [78.9]])

    def test_parse_int_field(self):
        field = fieldmetadata.FieldMetadata(
            'TCR',
            '1',
            'Integer',
            'Total reverse, strand coverage at this locus',
            'Blah',
            '0.3.0')
        self.assertEqual(field.name, "TCR")
        self.assertEqual(
            field.description,
            "Total reverse, strand coverage at this locus")
        self.assertEqual(field.source, "Blah")
        self.assertEqual(field.version, "0.3.0")
        self.assertEqual(field.extract_data("12", 1, None), [[12]])
        self.assertEqual(field.extract_data("12", 2, None), [[12], [12]])

    def test_parse_field_of_length_R(self):
        field = fieldmetadata.FieldMetadata(None, 'R', 'Integer', None)
        self.assertEqual(field.extract_data("12,13", 1, None), [[12, 13]])
        self.assertEqual(field.extract_data(
            "12,13,14", 2, None), [[12, 13], [12, 14]])
        with self.assertRaises(AssertionError):
            field.extract_data("12,78", 2, None)

    def test_parse_field_of_length_G(self):
        field = fieldmetadata.FieldMetadata(None, 'G', 'Integer', None)
        self.assertEqual(field.extract_data("12,13", None, 2), [[12], [13]])
        self.assertEqual(
            field.extract_data(
                "12,13,14", None, 3), [
                [12], [13], [14]])
        with self.assertRaises(AssertionError):
            field.extract_data("12,78", None, 1)

    def test_parse_field_of_unknown_number(self):
        field = fieldmetadata.FieldMetadata(None, '.', 'Integer', None)
        self.assertEqual(field.extract_data("12,13", 1, None), [[12, 13]])

    def test_parse_flag_field(self):
        field = fieldmetadata.FieldMetadata(None, '.', 'Flag', None)
        self.assertEqual(True, field('1'))
        self.assertEqual(True, field('YES'))
        self.assertEqual(True, field('TRUE'))
        self.assertEqual(False, field('FALSE'))
        with self.assertRaises(weCallException):
            field("blah")


class TestMakeSplitSampleAltFunc(unittest.TestCase):
    def test_split_genotype_likelihood_with_correct_number_of_genotypes_haploid(self):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual([[1.0, 2.0]], split_func(
            [1.0, 2.0], 1, GenotypeCall("0")))

    def test_split_genotype_likelihood_with_correct_number_of_genotypes_diploid(self):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual([[1.0, 2.0, 3.0]], split_func(
            [1.0, 2.0, 3.0], 1, GenotypeCall("0/1")))

    @testfixtures.log_capture()
    def test_split_genotype_likelihood_with_missing_genotype_likelihood_haploid(self, log):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual(
            [[None, None], [None, None]],
            split_func([1.0, 2.0], 2, GenotypeCall("0"))
        )
        log.check(
            ('wecall.vcfutils.fieldmetadata', 'WARNING',
             "Incorrect number of values 'G' cardinality, expected 3, got 2"),
        )

    @testfixtures.log_capture()
    def test_split_genotype_likelihood_with_missing_genotype_likelihood_diploid(self, log):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual(
            [[None, None, None], [None, None, None]],
            split_func([1.0, 2.0, 3.0], 2, GenotypeCall("0/1"))
        )
        log.check(
            ('wecall.vcfutils.fieldmetadata', 'WARNING',
             "Incorrect number of values 'G' cardinality, expected 6, got 3"),
        )

    def test_split_genotype_likelihood_with_correct_number_of_genotypes_haploid_multi_allelic(self):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual(
            [[1.0, 2.0], [1.0, 3.0]],
            split_func([1.0, 2.0, 3.0], 2, GenotypeCall("0"))
        )

    def test_split_genotype_likelihood_with_correct_number_of_genotypes_diploid_multi_allelic(self):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual(
            [[1.0, 2.0, 3.0], [1.0, 4.0, 6.0]],
            split_func([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 2, GenotypeCall("0/1"))
        )

    @testfixtures.log_capture()
    def test_split_genotype_likelihood_warns_for_no_genotype(self, log):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual(
            [[1.0, 2.0], [1.0, 2.0]],
            split_func([1.0, 2.0], 2, None)
        )
        log.check(
            ('wecall.vcfutils.fieldmetadata', 'WARNING',
             "Unknown ploidy when parsing genotype likelihood"),
        )

    @testfixtures.log_capture()
    def test_split_genotype_likelihood_warns_for_non_haploid_diploid(self, log):
        split_func = make_split_sample_alt_func("G", lambda x: x)
        self.assertEqual(
            [[1.0, 2.0], [1.0, 2.0]],
            split_func([1.0, 2.0], 2, GenotypeCall("0/1/2"))
        )
        log.check(
            ('wecall.vcfutils.fieldmetadata', 'WARNING',
             "Unable to handle ploidy other than haploid or diploid."),
        )
