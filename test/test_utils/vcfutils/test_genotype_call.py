# All content Copyright (C) 2018 Genomics plc
import unittest
from wecall.vcfutils.genotype_call import GenotypeCall, merge_genotype_calls
from wecall.vcfutils.sample_data import SampleData


class TestGenotypeCall(unittest.TestCase):
    def test_should_parse_alleles_into_sorted_list(self):
        self.assertEqual([None, 1, 2], GenotypeCall("2/./1").alleles)

    def test_should_report_when_genotype_is_unknown(self):
        genotype_call1 = GenotypeCall("././.")
        genotype_call2 = GenotypeCall("./0/.")
        self.assertTrue(genotype_call1.is_unknown())
        self.assertFalse(genotype_call2.is_unknown())

    def test_should_compare_equal_if_alleles_match_setwise_if_both_unphased(self):
        genotype_call1 = GenotypeCall("0/1")
        genotype_call2 = GenotypeCall("1/0")
        self.assertEqual(genotype_call1, genotype_call2)

    def test_should_compare_not_equal_if_alleles_match_setwise_if_both_phased(self):
        heterozygous_genotype_call1 = GenotypeCall("0|1")
        heterozygous_genotype_call2 = GenotypeCall("1|0")
        self.assertNotEqual(
            heterozygous_genotype_call1,
            heterozygous_genotype_call2)

    def test_should_compare_equal_if_alleles_match_setwise_if_one_phased_and_one_unphased(self):
        heterozygous_genotype_call_unphased = GenotypeCall("0/1")
        heterozygous_genotype_call_phased = GenotypeCall("1|0")
        self.assertEqual(
            heterozygous_genotype_call_unphased,
            heterozygous_genotype_call_phased)

    def test_should_compare_equal_if_all_alleles_match_regardless_of_phasing(self):
        genotype_call_with_phase = GenotypeCall("1|1")
        genotype_call_without_phase = GenotypeCall("1/1")
        self.assertEqual(genotype_call_with_phase, genotype_call_without_phase)

    def test_should_not_compare_equal_diploid_and_haploid_ref_call(self):
        diploid_call = GenotypeCall("0/0")
        haploid_call = GenotypeCall("0")
        self.assertNotEqual(diploid_call, haploid_call)

    def test_should_not_compare_equal_diploid_and_haploid_hom_call(self):
        diploid_call = GenotypeCall("1/1")
        haploid_call = GenotypeCall("1")
        self.assertNotEqual(diploid_call, haploid_call)


class TestGenotypeCallClassifications(unittest.TestCase):
    def test_should_mark_following_as_heterozygous(self):
        self.assertTrue(GenotypeCall("0/1").is_heterozygous())
        self.assertTrue(GenotypeCall("1/0").is_heterozygous())
        self.assertTrue(GenotypeCall("1/.").is_heterozygous())
        self.assertTrue(GenotypeCall("./1").is_heterozygous())
        self.assertTrue(GenotypeCall("0|1").is_heterozygous())
        self.assertTrue(GenotypeCall("1|0").is_heterozygous())
        self.assertTrue(GenotypeCall("1|.").is_heterozygous())
        self.assertTrue(GenotypeCall(".|1").is_heterozygous())
        self.assertTrue(GenotypeCall("1|2").is_heterozygous())

    def test_should_mark_following_as_not_heterozygous(self):
        self.assertFalse(GenotypeCall("./.").is_heterozygous())
        self.assertFalse(GenotypeCall(".|.").is_heterozygous())
        self.assertFalse(GenotypeCall("1/1").is_heterozygous())
        self.assertFalse(GenotypeCall("1|1").is_heterozygous())
        self.assertFalse(GenotypeCall("2/2").is_heterozygous())
        self.assertFalse(GenotypeCall("2|2").is_heterozygous())

    def test_should_mark_as_homozygous_alt(self):
        self.assertTrue(GenotypeCall("1/1").is_homozygous_alt())
        self.assertTrue(GenotypeCall("1|1").is_homozygous_alt())
        self.assertTrue(GenotypeCall("2/2").is_homozygous_alt())
        self.assertTrue(GenotypeCall("2|2").is_homozygous_alt())

    def test_should_not_mark_as_homozygous_alt(self):
        self.assertFalse(GenotypeCall("./.").is_homozygous_alt())
        self.assertFalse(GenotypeCall(".|.").is_homozygous_alt())
        self.assertFalse(GenotypeCall("0/1").is_homozygous_alt())
        self.assertFalse(GenotypeCall("1/0").is_homozygous_alt())
        self.assertFalse(GenotypeCall("1/.").is_homozygous_alt())
        self.assertFalse(GenotypeCall("./1").is_homozygous_alt())
        self.assertFalse(GenotypeCall("0|1").is_homozygous_alt())
        self.assertFalse(GenotypeCall("1|0").is_homozygous_alt())
        self.assertFalse(GenotypeCall("1|.").is_homozygous_alt())
        self.assertFalse(GenotypeCall(".|1").is_homozygous_alt())
        self.assertFalse(GenotypeCall("1|2").is_homozygous_alt())

    def test_should_mark_following_as_called(self):
        self.assertTrue(GenotypeCall("0/1").is_called())
        self.assertTrue(GenotypeCall("0|1").is_called())
        self.assertTrue(GenotypeCall("./1").is_called())
        self.assertTrue(GenotypeCall(".|1").is_called())
        self.assertTrue(GenotypeCall("0/2").is_called())
        self.assertTrue(GenotypeCall("1/2").is_called())
        self.assertTrue(GenotypeCall("././1").is_called())
        self.assertTrue(GenotypeCall("0/0/1").is_called())

    def test_should_not_mark_following_as_called(self):
        self.assertFalse(GenotypeCall("./.").is_called())
        self.assertFalse(GenotypeCall("./0").is_called())
        self.assertFalse(GenotypeCall("0/.").is_called())
        self.assertFalse(GenotypeCall("0|.").is_called())
        self.assertFalse(GenotypeCall(".|0").is_called())
        self.assertFalse(GenotypeCall(".").is_called())
        self.assertFalse(GenotypeCall("0").is_called())
        self.assertFalse(GenotypeCall("././.").is_called())
        self.assertFalse(GenotypeCall(".|0|.").is_called())


class TestMergeGenotypeCall(unittest.TestCase):
    # Arithmetic for GenotypeCall.
    # "0/0" + "0/0" = "0/0"
    # "0/0" + other = other
    # "1/1" + "1/1" = "1/1"
    # "1/1" + A => error for A in {0/1, 0|1, 1|0}
    # "0/1" + "0/1" => "1/1" (two calls made at different points.)
    # "0/1" + "1|0" = "0/1" + "0|1" = "1/1" (??)
    # "0|1" + "1|0" = "1/1"
    # "0|1" + "0|1" = "0|1"
    # "1|0" + "1|0" = "1|0"

    def test_should_get_homozygous_alt_if_combining_two_homozygous_alt_genotypes(self):
        genotype_call_1 = GenotypeCall("1/1")
        genotype_call_2 = GenotypeCall("1/1")
        self.assertEqual(
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2),
            GenotypeCall("1/1"))

    def test_should_get_homozygous_ref_if_combining_two_homozygous_ref_genotypes(self):
        genotype_call_1 = GenotypeCall("0/0")
        genotype_call_2 = GenotypeCall("0/0")
        self.assertEqual(
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2),
            GenotypeCall("0/0"))

    def test_should_raise_exception_one_is_homozygous_alt_and_other_heterozygous_genotypes(self):
        genotype_call_1 = GenotypeCall("0/1")
        genotype_call_2 = GenotypeCall("1/1")
        self.assertRaises(
            Exception,
            merge_genotype_calls,
            genotype_call_1,
            genotype_call_2)

    def test_should_get_heterozygous_if_one_is_homozygous_ref_and_other_is_heterozygous(self):
        genotype_call_1 = GenotypeCall("0/0")
        genotype_call_2 = GenotypeCall("0/1")
        self.assertEqual(
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2),
            genotype_call_2)

    def test_should_get_homozygous_alt_if_one_is_homozyzgous_ref_and_other_is_homozygous_alt(self):
        genotype_call_1 = GenotypeCall("0/0")
        genotype_call_2 = GenotypeCall("1/1")
        self.assertEqual(
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2),
            genotype_call_2)

    def test_should_combine_two_unphased_heterozygous_genotypes_to_homozygous_alt(self):
        genotype_call_1 = GenotypeCall("0/1")
        genotype_call_2 = GenotypeCall("0/1")
        self.assertEqual(
            GenotypeCall("1/1"),
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2))

    def test_should_combine_two_heterozygous_genotypes_to_homozygous_alt_if_one_is_not_phased(self):
        genotype_call_1 = GenotypeCall("1|0")
        genotype_call_2 = GenotypeCall("0/1")
        self.assertEqual(
            GenotypeCall("1/1"),
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2))

        genotype_call_1 = GenotypeCall("0|1")
        genotype_call_2 = GenotypeCall("0/1")
        self.assertEqual(
            GenotypeCall("1/1"),
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2))

    def test_should_return_phased_heterozygous_genotype_when_merging_two_phased_identical_heterozygous_genotypes(self):
        genotype_call_1 = GenotypeCall("1|0")
        genotype_call_2 = GenotypeCall("1|0")
        self.assertEqual(
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2),
            GenotypeCall("1|0"))

        genotype_call_1 = GenotypeCall("0|1")
        genotype_call_2 = GenotypeCall("0|1")
        self.assertEqual(
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2),
            GenotypeCall("0|1"))

    def test_should_combine_two_opposite_heterozygous_phased_genotypes(self):
        genotype_call_1 = GenotypeCall("1|0")
        genotype_call_2 = GenotypeCall("0|1")
        self.assertEqual(
            GenotypeCall("1|1"),
            merge_genotype_calls(
                genotype_call_1,
                genotype_call_2))

    def test_should_raise_exception_any_genotype_is_not_diploid(self):
        genotype_call_1 = GenotypeCall("1")
        genotype_call_2 = GenotypeCall("0")
        self.assertRaises(
            Exception,
            merge_genotype_calls,
            genotype_call_1,
            genotype_call_2)


class TestMergeSampleDataGenotypesCalls(unittest.TestCase):

    def test_should_fail_if_sample_data_objects_have_different_sample(self):
        sample_data1 = SampleData(['GT'], ['sample_name_1'])
        sample_data1.add_sample_data(
            'sample_name_1', 'GT', GenotypeCall('0/0'))
        sample_data2 = SampleData(['GT'], ['sample_name_2'])
        sample_data2.add_sample_data(
            'sample_name_2', 'GT', GenotypeCall('0/0'))

        self.assertRaises(
            Exception,
            sample_data1.merge_genotype_calls,
            sample_data2.genotypes())

    def test_should_merge_genotype_call_object_in_sample_data(self):
        sample_data1 = SampleData(['GT'], ['sample_name'])
        sample_data1.add_sample_data('sample_name', 'GT', GenotypeCall('0/1'))
        sample_data2 = SampleData(['GT'], ['sample_name'])
        sample_data2.add_sample_data('sample_name', 'GT', GenotypeCall('0/1'))

        sample_data1.merge_genotype_calls(sample_data2.genotypes())

        self.assertEqual(
            sample_data1.get_field("sample_name", "GT"),
            GenotypeCall("1/1")
        )


class TestGenotypeCallAlleleCounts(unittest.TestCase):

    def test_homozygous_unphased_genotypes(self):
        self.assertEqual((1,), GenotypeCall('0').normalized_allele_count)
        self.assertEqual((1,), GenotypeCall('0/0').normalized_allele_count)
        self.assertEqual((1,), GenotypeCall('0/0/0').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('1').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('1/1').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('1/1/1').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('2').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('2/2').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('2/2/2').normalized_allele_count)

    def test_homozygous_phased_genotypes(self):
        self.assertEqual((1,), GenotypeCall('0').normalized_allele_count)
        self.assertEqual((1,), GenotypeCall('0|0').normalized_allele_count)
        self.assertEqual((1,), GenotypeCall('0|0|0').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('1').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('1|1').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('1|1|1').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('2').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('2|2').normalized_allele_count)
        self.assertEqual((0, 1), GenotypeCall('2|2|2').normalized_allele_count)

    def test_binary_heterozygous_unphased_genotypes(self):
        self.assertEqual((1, 1), GenotypeCall('0/1').normalized_allele_count)
        self.assertEqual((1, 1), GenotypeCall('0/2').normalized_allele_count)
        self.assertEqual((0, 1, 1), GenotypeCall(
            '1/2').normalized_allele_count)

    def test_unknown_genotypes_allele_count(self):
        self.assertEqual((1, ), GenotypeCall('.').normalized_allele_count)
        self.assertEqual((1, ), GenotypeCall('./.').normalized_allele_count)
        self.assertEqual((1, 1), GenotypeCall('./1').normalized_allele_count)

    def test_binary_heterozygous_phased_genotypes(self):
        self.assertEqual((1, 1), GenotypeCall('0|1').normalized_allele_count)
        self.assertEqual((1, 1), GenotypeCall('1|0').normalized_allele_count)
        self.assertEqual((1, 1), GenotypeCall('0|2').normalized_allele_count)
        self.assertEqual((1, 1), GenotypeCall('2|0').normalized_allele_count)
        self.assertEqual((0, 1, 1), GenotypeCall(
            '1|2').normalized_allele_count)
        self.assertEqual((0, 1, 1), GenotypeCall(
            '2|1').normalized_allele_count)
