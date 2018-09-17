# All content Copyright (C) 2018 Genomics plc
from collections import OrderedDict
import copy
from wecall.bamutils.sample_bank import SampleBank
from wecall.bamutils.sequence import Sequence
from wecall.common.exceptions import EchidnaException
from wecall.genomics.variant import Variant
from wecall.vcfutils.genotype_call import GenotypeCall
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.variant_caller_builder import VariantCallerBuilderFromSampleBank

DEFAULT_SAMPLE_NAME = 'TEST_SAMPLE'


class AsciiWecallRunnerTest(BaseTest):
    def calls_variants(
            self,
            ref,
            sequence_list,
            expected_ascii_haplotypes=None,
            expected_variant_stubs=None,
            n_fwd=None,
            n_rev=None,
            config_dict=None, ):
        sample_bank = self.__build_default_sample_bank(
            ref, sequence_list, n_fwd, n_rev)
        variant_callset = self.__run_wecall(sample_bank, config_dict)
        expected_variants = self._get_expected_variants(
            expected_ascii_haplotypes, expected_variant_stubs, sample_bank)

        self.assertEqual(variant_callset.get_variants(), expected_variants)

    def calls_variants_with_genotype(
            self,
            ref,
            sequence_list,
            expected_haplotypes=None,
            expected_variants_with_genotypes=None,
            config_dict=None):
        self.__validate_expected_calls(
            expected_haplotypes,
            expected_variants_with_genotypes)
        sample_bank = self.__build_default_sample_bank(ref, sequence_list)
        variant_callset = self.__run_wecall(sample_bank, config_dict)
        wecall_calls = variant_callset.get_variants_with_genotypes()

        if expected_variants_with_genotypes is None:
            expected_calls_for_default_sample = {
                sample_bank.sample_names[0]: expected_haplotypes}
            expected_calls = self.__get_expected_calls_from_sample_ascii_haplotypes(
                expected_calls_for_default_sample, sample_bank.reference)
        else:
            expected_calls = {}
            sample_name = sample_bank.sample_names[0]
            for expected_stub in expected_variants_with_genotypes:
                variant = self._variant_from_stub(
                    sample_bank.reference.chrom, expected_stub)
                expected_calls[variant] = {
                    sample_name: GenotypeCall(
                        expected_stub[3])}

        self.assertDictEqual(expected_calls, wecall_calls)

    def calls_variants_with_info_annotation(
            self,
            ref,
            sequence_list,
            expected_variants_and_info_annotation,
            config_dict=None,
    ):
        sample_bank = self.__build_default_sample_bank(ref, sequence_list)
        variant_callset = self.__run_wecall(sample_bank, config_dict)
        expected_calls = {}

        for expected_stub in expected_variants_and_info_annotation:
            variant = self._variant_from_stub(
                sample_bank.reference.chrom, expected_stub)
            expected_calls[variant] = {}

            for sample_data_key, value in expected_stub[3].items():
                expected_calls[variant][sample_data_key] = value

        actual_calls = variant_callset.get_variant_records()

        # assert the expected and actual variants called are the same
        self.assertEqual(set(expected_calls.keys()), set(actual_calls.keys()))

        for variant in expected_calls:
            expected_sample_data = expected_calls[variant]
            actual_sample_data = {}

            for field in expected_sample_data:
                value = actual_calls[variant].info[field]
                actual_sample_data[field] = value

            self.assertEqual(expected_sample_data, actual_sample_data)

    def __run_wecall(self, sample_bank, config_dict):
        vc_builder = VariantCallerBuilderFromSampleBank(
            sample_bank, self.work_dir)
        if config_dict is not None:
            for key, value in config_dict.items():
                vc_builder.configuration[key] = value
        vc_wrapper = vc_builder.build()
        vc_callset = vc_wrapper.run().get_variant_callset(self)
        return vc_callset

    def calls_variants_from_samples(self,
                                    ref,
                                    sample_seqs,
                                    expected_haplotypes=None,
                                    expected_call_stubs=None,
                                    config_dict=None):
        """
        :param expected_haplotypes: dictionary: {sample_name : list of two ascii sequences expressing the genotype}
        :param expected_call_stubs: dictionary: {variant_stub: dictionary {sample_name: str(genotype)} }
        """
        self.__validate_expected_calls(
            expected_haplotypes, expected_call_stubs)
        sample_bank = SampleBank(ref)

        for sample_name, sequence_list in sample_seqs.items():
            sample_bank.add_sample_with_seqs_and_quals(
                sample_name, sequence_list)

        variant_callset = self.__run_wecall(sample_bank, config_dict)
        wecall_calls = variant_callset.get_variants_with_genotypes()

        if expected_call_stubs is None:
            self.__filter_none_genotypes(wecall_calls)
            expected_calls = self.__get_expected_calls_from_sample_ascii_haplotypes(
                expected_haplotypes, sample_bank.reference)
        else:
            expected_calls = {}
            for variant_stub, genotypes in expected_call_stubs.items():
                variant = self._variant_from_stub(
                    sample_bank.reference.chrom, variant_stub)
                expected_calls[variant] = OrderedDict()
                for sample_name, genotype in genotypes.items():
                    expected_calls[variant][sample_name] = GenotypeCall(
                        genotype)

        self.maxDiff = None  # print the whole message if the following assertion fails
        self.assertDictEqual(expected_calls, wecall_calls)

    def calls_variants_with_sample_data_and_filters(
            self,
            ref,
            sequence_list,
            expected_variants_and_sample_data,
            config_dict=None):

        sample_bank = self.__build_default_sample_bank(ref, sequence_list)
        variant_callset = self.__run_wecall(sample_bank, config_dict)
        expected_calls = {}

        for expected_stub in expected_variants_and_sample_data:
            variant = self._variant_from_stub(
                sample_bank.reference.chrom, expected_stub)
            expected_calls[variant] = {}

            for sample_data_key, value in expected_stub[3].items():
                expected_calls[variant][sample_data_key] = value

            expected_calls[variant]["FILTERS"] = set()

            for filterValue in expected_stub[4]:
                if filterValue != "PASS":
                    expected_calls[variant]["FILTERS"].add(filterValue)

        sample_name = sample_bank.sample_names[0]
        actual_calls = variant_callset.get_variant_records()

        # assert the expected and actual variants called are the same
        self.assertEqual(set(expected_calls.keys()), set(actual_calls.keys()))

        for variant in expected_calls:
            expected_sample_data = expected_calls[variant]
            actual_sample_data = {}

            for field in expected_sample_data:
                if field != "FILTERS":
                    value = actual_calls[variant].sample_info.get_genotype_data(sample_name)[
                        field]
                    actual_sample_data[field] = value

            actual_sample_data["FILTERS"] = actual_calls[variant].filters

            self.assertEqual(expected_sample_data, actual_sample_data)

    def calls_variants_with_sample_data(self,
                                        ref,
                                        sequence_list,
                                        expected_variants_and_sample_data,
                                        config_dict=None):
        sample_bank = self.__build_default_sample_bank(ref, sequence_list)
        variant_callset = self.__run_wecall(sample_bank, config_dict)
        expected_calls = {}

        for expected_stub in expected_variants_and_sample_data:
            variant = self._variant_from_stub(
                sample_bank.reference.chrom, expected_stub)
            expected_calls[variant] = {}

            for sample_data_key, value in expected_stub[3].items():
                expected_calls[variant][sample_data_key] = value

        sample_name = sample_bank.sample_names[0]
        actual_calls = variant_callset.get_variant_records()

        # assert the expected and actual variants called are the same
        self.assertEqual(set(expected_calls.keys()), set(actual_calls.keys()))

        for variant in expected_calls:
            expected_sample_data = expected_calls[variant]
            actual_sample_data = {}

            for field in expected_sample_data:
                value = actual_calls[variant].sample_info.get_genotype_data(sample_name)[
                    field]
                actual_sample_data[field] = value

            self.assertEqual(expected_sample_data, actual_sample_data)

    def __build_default_sample_bank(
            self,
            ref,
            sequence_list,
            n_fwd=None,
            n_rev=None):
        sample_bank = SampleBank(ref)
        sample_bank.add_sample_with_seqs_and_quals(
            DEFAULT_SAMPLE_NAME, sequence_list, n_fwd, n_rev)

        return sample_bank

    @staticmethod
    def __validate_expected_calls(expected_ascii, expected_stubs):
        if expected_ascii is None and expected_stubs is None:
            raise EchidnaException(
                "Expected variants have to be provided either in the ascii or variant stub format."
            )

    @staticmethod
    def _get_expected_variants(
            ascii_haplotypes,
            expected_variant_stubs,
            sample_bank):
        if ascii_haplotypes is None and expected_variant_stubs is None:
            return sample_bank.variants
        elif ascii_haplotypes is None:
            return {
                AsciiWecallRunnerTest._variant_from_stub(
                    sample_bank.reference.chrom,
                    stub) for stub in expected_variant_stubs}
        else:
            return set(
                AsciiWecallRunnerTest.__get_expected_calls_from_haplotypes(
                    ascii_haplotypes,
                    sample_bank.reference).keys())

    @staticmethod
    def __get_expected_calls_from_sample_ascii_haplotypes(
            ascii_haplotypes, reference):
        calls_per_variant = {}
        for sample_name, ascii_strings in ascii_haplotypes.items():
            calls_for_sample = AsciiWecallRunnerTest.__get_expected_calls_from_haplotypes(
                ascii_strings, reference)
            for variant, genotype in calls_for_sample.items():
                if variant in calls_per_variant and sample_name in calls_per_variant[variant]:
                    raise EchidnaException(
                        "Cannot supply multiple genotypes for "
                        "sample_name {} and variant {}.".format(
                            sample_name, variant))
                if variant not in calls_per_variant:
                    # ordered dict only to comply with what the actual calls
                    # look like
                    calls_per_variant[variant] = OrderedDict()

                calls_per_variant[variant][sample_name] = genotype

        return calls_per_variant

    @staticmethod
    def __get_expected_calls_from_haplotypes(ascii_strings, reference):
        if len(ascii_strings) != 2:
            raise EchidnaException(
                "Expected calls have to be defined as a diploid.")
        if not all(len(str) == reference.length_with_deletions()
                   for str in ascii_strings):
            raise EchidnaException(
                "Ascii haplotypes have to be of the same length as the reference")

        vars_from_hap1 = Sequence(reference, ascii_strings[0]).variants
        vars_from_hap2 = Sequence(reference, ascii_strings[1]).variants

        calls = {}
        for var in vars_from_hap1.intersection(vars_from_hap2):
            calls[var] = GenotypeCall("1/1")
        for var in vars_from_hap1.symmetric_difference(vars_from_hap2):
            calls[var] = GenotypeCall("0/1")

        return calls

    @staticmethod
    def _variant_from_stub(chrom, stub):
        pos = stub[0]
        ref = stub[1]
        alt = stub[2]
        return Variant(chrom, pos, ref, alt)

    @staticmethod
    def __filter_none_genotypes(calls):
        calls_copy = copy.deepcopy(calls)
        for variant, genotypes in calls_copy.items():
            for sample_name, genotype in genotypes.items():
                if genotype == GenotypeCall("./."):
                    calls[variant].popitem(sample_name)
