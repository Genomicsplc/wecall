# All content Copyright (C) 2018 Genomics plc
import os
from wecall.genomics.variant import Variant
from wecall.utils.interval import ChromInterval
from wecall.vcfutils.genotype_call import GenotypeCall
from wecall.vcfutils.parser import VCFReaderContextManager

ref_alt = "<NON_REF>"


class VCFExpectation(object):

    def __init__(self, test_case, path):
        self.__test_case = test_case
        self.__path = path
        self.__test_case.assertTrue(os.path.exists(self.__path))
        with VCFReaderContextManager(self.__path) as vcf_reader:
            self.__schema = vcf_reader.read_header()
            self.__records = list(vcf_reader.read_records())

        # ensure genotype likelihoods are within range throughout
        for record in self.__records:
            for sample_name in record.sample_info.get_sample_names():
                try:
                    for GL_value in record.sample_info.get_field(
                            sample_name, 'GL'):
                        self.__test_case.assertTrue(GL_value <= 0.0)
                except KeyError:
                    pass

    def __eq__(self, other):
        return self.__records == other.__records

    def __ne__(self, other):
        return not self.__eq__(other)

    def record_count(self, expected):
        self.__test_case.assertEqual(expected, len(self.__records))
        return self

    def has_info_meta_data(self, key):
        self.__test_case.assertIn(
            key, {key for key, value in self.__schema.iter_info_data()})
        return VCFInfoMetadataExpectation(
            self.__test_case, self.__schema.get_info_data(key))

    def has_filter(self, key):
        self.__test_case.assertIn(
            key, {key for key, value in self.__schema.iter_filters()})
        return VCFFilterMetadataExpectation(
            self.__test_case, self.__schema.get_filter(key))

    def has_reference_calls_for_region(self, chrom, start, end):
        return self.has_reference_calls(ChromInterval(chrom, start, end))

    def has_reference_calls(self, chrom_interval):
        records = {record.variant: record for record in self.__records if record.variant.overlap(
            chrom_interval)}
        # check all are ref calls.
        current_position = chrom_interval.start
        for variant in sorted(records.keys()):
            self.__test_case.assertEqual(variant.pos_from, current_position)
            self.__test_case.assertEqual(variant.alt, ref_alt)
            record = records[variant]
            info_expectation = VCFRecordExpectation(
                self.__test_case, record).with_info()

            info_expectation.with_field("BEG", [current_position + 1])
            current_end = record.info["END"][0]
            info_expectation.with_field(
                "LEN", [current_end - current_position])

            # Due to mix of closed intervals and 1 --> 0 indexing changes this
            # is true!
            current_position = current_end
        self.__test_case.assertEqual(chrom_interval.end, current_position)
        return self

    def has_record(self, chrom, pos, ref, alt):
        return self.has_record_for_variant(Variant(chrom, pos, ref, alt))

    def has_record_for_variant(self, variant):
        records = {record.variant: record for record in self.__records}
        self.__test_case.assertIn(variant, records.keys())
        return VCFRecordExpectation(self.__test_case, records[variant])

    def missing_record_for_variant(self, variant):
        self.__test_case.assertNotIn(
            variant, {record.variant for record in self.__records})
        return self

    def has_record_for_variants(self, *variants):
        records = {record.variant: record for record in self.__records}
        records_to_check = list()
        for variant in variants:
            self.__test_case.assertIn(variant, records)
            records_to_check.append(records[variant])
        return VCFRecordListExpectation(
            self.__test_case, [
                VCFRecordExpectation(
                    self.__test_case, record) for record in records_to_check])

    def with_samples(self, sample_names):
        self.__test_case.assertEqual(sample_names, self.__schema.samples)
        return self


class VCFRecordListExpectation(object):
    def __init__(self, test_case, records):
        self.__test_case = test_case
        self.__record_expectations = records

    def with_sample(self, sample_name):
        sample_expectations = [record_expecation.with_sample(
            sample_name) for record_expecation in self.__record_expectations]
        return VCFSampleDataListExpectation(
            self.__test_case, sample_expectations)


class VCFRecordExpectation(object):
    def __init__(self, test_case, record):
        self.__test_case = test_case
        self.__record = record

    def with_sample(self, sample_name):
        self.__test_case.assertTrue(
            self.__record.sample_info.has_sample(sample_name))
        return VCFSampleDataExpectation(
            self.__test_case,
            self.__record.sample_info,
            sample_name)

    def with_from_multi_allelic_record(self, expected=True):
        self.__test_case.assertEqual(expected, self.__record.from_multi_alt)
        return self

    def with_sample_info(self, expected_sample_info):
        self.__test_case.assertEqual(
            expected_sample_info,
            self.__record.sample_info)
        return self

    def with_info_data(self, expected_info):
        self.__test_case.assertEqual(expected_info, self.__record.info)
        return self

    def with_quality(self, expected_quality):
        self.__test_case.assertEqual(expected_quality, self.__record.quality)
        return self

    def with_no_filters(self):
        return self.with_filters(set())

    def with_filters(self, expected_filters):
        self.__test_case.assertEqual(expected_filters, self.__record.filters)
        return self

    def with_ids(self, expected_ids):
        self.__test_case.assertEqual(expected_ids, self.__record.ids)
        return self

    def with_info(self):
        return VCFInfoFieldExpectation(self.__test_case, self.__record.info)


class VCFInfoFieldExpectation(object):
    def __init__(self, test_case, info_data):
        self.__test_case = test_case
        self.__info_data = info_data

    def with_field(self, key, values):
        self.__test_case.assertEqual(values, self.__info_data[key])
        return self

    def has_keys(self, *expected_keys):
        self.__test_case.assertEqual(
            set(expected_keys), set(
                self.__info_data.keys()))
        return self


class VCFFilterMetadataExpectation(object):

    def __init__(self, test_case, filter_meta_data):
        self.__test_case = test_case
        self.__filter_meta_data = filter_meta_data

    def has_description(self, expected_description):
        self.__test_case.assertEqual(
            expected_description,
            self.__filter_meta_data.description)
        return self


class VCFInfoMetadataExpectation(object):

    def __init__(self, test_case, info_meta_data):
        self.__test_case = test_case
        self.__info_meta_data = info_meta_data

    def has_description(self, expected_description):
        self.__test_case.assertEqual(
            expected_description,
            self.__info_meta_data.description)
        return self

    def has_data_type(self, expected_type):
        self.__test_case.assertEqual(
            expected_type, self.__info_meta_data.data_type)
        return self


class VCFSampleDataListExpectation(object):
    def __init__(self, test_case, sample_data_expectations):
        self.__test_case = test_case
        self.__sample_data_expectations = sample_data_expectations

    def has_phased_genotypes(self, *genotype_strings):
        actual_genotypes = tuple(str(sample_data_expectation.get_genotype(
        )) for sample_data_expectation in self.__sample_data_expectations)
        actual_reverse_genotypes = tuple(
            str(sample_data_expectation.get_genotype())[::-1]
            for sample_data_expectation in self.__sample_data_expectations
        )

        expected_genotypes = tuple(genotype_strings)
        expected_reverse_genotypes = tuple(
            genotype_string[::-1] for genotype_string in genotype_strings
        )
        assert(len(actual_genotypes) == len(expected_genotypes))
        self.__test_case.assertEqual(
            {expected_genotypes, expected_reverse_genotypes},
            {actual_genotypes, actual_reverse_genotypes}
        )
        return self

    def has_exact_phased_genotypes(self, *genotype_strings):
        actual_genotypes = tuple(str(sample_data_expectation.get_genotype(
        )) for sample_data_expectation in self.__sample_data_expectations)

        expected_genotypes = tuple(genotype_strings)
        assert(len(actual_genotypes) == len(expected_genotypes))
        self.__test_case.assertEqual(expected_genotypes, actual_genotypes)
        return self

    def has_phase_set_id(self, expected):
        for sample_data_expectation in self.__sample_data_expectations:
            sample_data_expectation.has_phase_set_id(expected)
        return self

    def has_phase_set_quality(self, expected):
        for sample_data_expectation in self.__sample_data_expectations:
            sample_data_expectation.has_phase_set_quality(expected)
        return self


class VCFSampleDataExpectation(object):

    phred_likelihood_key = 'PL'
    read_depth_key = 'DP'
    min_read_depth_key = 'MIN_DP'
    allelic_depth_key = 'AD'
    phase_set_key = 'PS'
    phase_set_quality = 'PQ'
    variant_allelic_frequency_key = 'VAF'

    def __init__(self, test_case, sample_data, sample_name):
        self.__test_case = test_case
        self.__sample_name = sample_name
        self.__sample_data = sample_data

    def get_genotype(self):
        genotypes = self.__sample_data.genotypes()
        self.__test_case.assertIn(self.__sample_name, genotypes)
        return genotypes[self.__sample_name]

    def has_phase_set_id(self, expected_id):
        self.__test_case.assertTrue(
            self.__sample_data.has_genotype_key(
                self.phase_set_key))
        self.__test_case.assertEqual(
            [expected_id], self.__sample_data.get_field(
                self.__sample_name, self.phase_set_key))
        return self

    def has_phase_set_quality(self, expected):
        self.__test_case.assertTrue(
            self.__sample_data.has_genotype_key(
                self.phase_set_quality))
        self.__test_case.assertEqual(
            [expected], self.__sample_data.get_field(
                self.__sample_name, self.phase_set_quality))
        return self

    def has_phased_genotype(self, genotype_string):
        actual_genotype_call = self.get_genotype()
        expected_genotype_call = GenotypeCall(genotype_string)

        self.__test_case.assertEqual(
            expected_genotype_call,
            actual_genotype_call)
        self.__test_case.assertEqual(
            expected_genotype_call.phased,
            actual_genotype_call.phased)
        return self

    def has_genotype(self, genotype_string):
        actual_genotype_call = self.get_genotype()
        expected_genotype_call = GenotypeCall(genotype_string)

        self.__test_case.assertEqual(
            expected_genotype_call,
            actual_genotype_call)
        return self

    def has_allelic_read_support(self, reference, *alts):
        self.__test_case.assertTrue(
            self.__sample_data.has_genotype_key(
                self.allelic_depth_key))
        ref_and_allelic_depths = self.__sample_data.get_field(
            self.__sample_name, self.allelic_depth_key)
        self.__test_case.assertEqual(
            len(ref_and_allelic_depths), 1 + len(alts))
        self.__test_case.assertEqual(ref_and_allelic_depths[0], reference)
        self.__test_case.assertEqual(ref_and_allelic_depths[1:], list(alts))
        return self

    def has_variant_allelic_frequency(self, *allelic_frequencies):
        self.__test_case.assertTrue(
            self.__sample_data.has_genotype_key(
                self.variant_allelic_frequency_key))
        actual_allele_frequencies = self.__sample_data.get_field(
            self.__sample_name,
            self.variant_allelic_frequency_key
        )
        self.__test_case.assertEqual(
            len(actual_allele_frequencies),
            len(allelic_frequencies))

        for actual_allele_frequency, expected_allele_frequency in zip(
                actual_allele_frequencies, allelic_frequencies):
            if expected_allele_frequency == ".":
                self.__test_case.assertEqual(
                    expected_allele_frequency, actual_allele_frequency)
            else:
                self.__test_case.assertAlmostEqual(
                    expected_allele_frequency, actual_allele_frequency, places=3)
        return self

    def has_read_depth(self, expected):
        self.__test_case.assertTrue(self.__sample_data.has_genotype_key(self.read_depth_key))
        read_depth = self.__sample_data.get_field(self.__sample_name, self.read_depth_key)
        self.__test_case.assertEqual(len(read_depth), 1)
        self.__test_case.assertEqual(read_depth[0], expected)
        return self

    def has_min_read_depth(self, expected):
        self.__test_case.assertTrue(self.__sample_data.has_genotype_key(self.min_read_depth_key))
        min_read_depth = self.__sample_data.get_field(self.__sample_name, self.min_read_depth_key)
        self.__test_case.assertEqual(len(min_read_depth), 1)
        self.__test_case.assertEqual(min_read_depth[0], expected)
        return self

    def has_RR_genotype_likelihood(self, reference_value):
        likelihoods = self.__sample_data.get_field(
            self.__sample_name, self.phred_likelihood_key)
        self.__test_case.assertEqual(len(likelihoods), 3)
        self.__test_case.assertEqual(likelihoods[0], reference_value)
        return self

    def has_RA_genotype_likelihood(self, reference_value):
        likelihoods = self.__sample_data.get_field(
            self.__sample_name, self.phred_likelihood_key)
        self.__test_case.assertEqual(len(likelihoods), 3)
        self.__test_case.assertEqual(likelihoods[1], reference_value)
        return self

    def has_AA_genotype_likelihood(self, reference_value):
        likelihoods = self.__sample_data.get_field(
            self.__sample_name, self.phred_likelihood_key)
        self.__test_case.assertEqual(len(likelihoods), 3)
        self.__test_case.assertEqual(likelihoods[2], reference_value)
        return self
