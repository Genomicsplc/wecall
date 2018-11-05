# All content Copyright (C) 2018 Genomics plc
from collections import OrderedDict

from wecall.common.exceptions import weCallException
from wecall.vcfutils.genotype_call import GenotypeCall, merge_genotype_calls
from wecall.vcfutils.stringutils import to_vcf_str


GENOTYPE_QUALITY_KEY = 'GQ'

GENOTYPE_KEY = 'GT'
GENOTYPE_LIKELIHOODS_KEY = 'GL'
GENOTYPE_PHRED_LIKELIHOODS_KEY = 'PL'
LIKELIHOOD_SCALING_FACTOR = OrderedDict(
    [(GENOTYPE_PHRED_LIKELIHOODS_KEY, -10.0), (GENOTYPE_LIKELIHOODS_KEY, 1.0)])
GENOTYPE_LIKELIHOODS_KEYS = list(LIKELIHOOD_SCALING_FACTOR.keys())

NUMBER_OF_READ = "NR"
READ_DEPTH = "DP"
READ_DEPTH_KEYS = [READ_DEPTH, NUMBER_OF_READ]

ALLELIC_DEPTH = "AD"
VARIANT_SUPPORT = "NV"
VARIANT_SUPPORT_MAP = {
    ALLELIC_DEPTH: lambda v: v[1:], VARIANT_SUPPORT: lambda v: v}


class SampleData(object):

    __slots__ = (
        '__sample_names',
        '__key_to_sample_values',
        '__merged_genotypes')

    def __init__(self, key_names, sample_names):
        self.__sample_names = sample_names
        self.__key_to_sample_values = OrderedDict()
        self.__merged_genotypes = False

        for key_name in key_names:
            default_value = []

            if key_name == GENOTYPE_KEY:
                default_value = GenotypeCall('./.')

            self.__key_to_sample_values[key_name] = [
                default_value for _ in range(len(self.__sample_names))]

    def __eq__(self, other):
        return (
            self.__sample_names == other.__sample_names and
            self.__key_to_sample_values == other.__key_to_sample_values
        )

    def __repr__(self):
        return "<SampleData: samples={!r}, keys={!r}>".format(
            self.__sample_names,
            list(self.__key_to_sample_values.items())
        )

    def add_sample_data(self, sample_name, key_name, sample_data_value):
        if key_name not in self.__key_to_sample_values:
            raise weCallException(
                "Missing key {} when adding sample data.".format(key_name))

        if sample_name not in self.__sample_names:
            raise weCallException(
                "Missing sample name {} supplied when adding sample data.".format(sample_name))

        if key_name == GENOTYPE_KEY and not isinstance(
                sample_data_value, GenotypeCall):
            raise weCallException("Genotype field must be a GenotypeCall.")

        self.__key_to_sample_values[key_name][self.__sample_names.index(
            sample_name)] = sample_data_value

    def get_genotype_data(self, sample_name):
        return SampleDataView(self, sample_name)

    def get_genotype_quality(self, sample_name):
        return self.get_field(sample_name, GENOTYPE_QUALITY_KEY)

    def has_read_depth_key(self):
        return any(self.has_genotype_key(key) for key in READ_DEPTH_KEYS)

    def get_read_depth(self, sample_name):
        for key in READ_DEPTH_KEYS:
            if self.has_genotype_key(key):
                return self.get_field(sample_name, key)
        raise weCallException(
            "Expected one of {} as the depth key.".format(READ_DEPTH_KEYS))

    def has_variant_support_key(self):
        return any(self.has_genotype_key(key)
                   for key in list(VARIANT_SUPPORT_MAP.keys()))

    def get_variant_support(self, sample_name):
        for key in list(VARIANT_SUPPORT_MAP.keys()):
            if self.has_genotype_key(key):
                return VARIANT_SUPPORT_MAP[key](
                    self.get_field(sample_name, key))
        raise weCallException(
            "Expected one of {} as the variant support key.".format(
                list(
                    VARIANT_SUPPORT_MAP.keys())))

    def has_no_likelihoods(self):
        if self.has_genotype_key(GENOTYPE_LIKELIHOODS_KEY) or self.has_genotype_key(
                GENOTYPE_PHRED_LIKELIHOODS_KEY):
            for sample_name in self.__sample_names:
                likelihoods = self.get_genotype_likelihoods(sample_name)
                if any([likeli is not None for likeli in likelihoods]):
                    return False
        return True

    def get_raw_genotype_likelihoods(self, sample_name):
        for key in list(LIKELIHOOD_SCALING_FACTOR.keys()):
            if self.has_genotype_key(key):
                try:
                    return self.get_field(sample_name, key)
                except KeyError:
                    pass
        return None

    def get_genotype_likelihoods(self, sample_name):
        def convert_likelihoods(likelihoods, factor):
            if likelihoods is None or likelihoods == '.':
                return likelihoods
            else:
                return [
                    None if value in {
                        None,
                        '.'} else value /
                    factor for value in likelihoods]

        for key in list(LIKELIHOOD_SCALING_FACTOR.keys()):
            if self.has_genotype_key(key):
                values = self.get_field(sample_name, key)
                return convert_likelihoods(
                    values, LIKELIHOOD_SCALING_FACTOR[key])
        raise weCallException(
            "Expected one of {} as the likelihood key.".format(
                list(
                    LIKELIHOOD_SCALING_FACTOR.keys())))

    def set_genotype_likelihoods(self, sample_name, likelihood_values):
        def convert_likelihoods(likelihoods, factor):
            if likelihoods is None or likelihoods == '.':
                return likelihoods
            else:
                return [
                    None if value in {
                        None,
                        '.'} else value *
                    factor for value in likelihoods]

        for key in list(LIKELIHOOD_SCALING_FACTOR.keys()):
            if self.has_genotype_key(key):
                converted_values = convert_likelihoods(
                    likelihood_values, LIKELIHOOD_SCALING_FACTOR[key])
                self.add_sample_data(sample_name, key, converted_values)
                return
        raise weCallException(
            "Expected one of {} as the likelihood key.".format(
                list(
                    LIKELIHOOD_SCALING_FACTOR.keys())))

    def get_field(self, sample_name, key):
        return self.__key_to_sample_values[key][self.__sample_names.index(
            sample_name)]

    def get_fields(self, sample_name):
        index = self.__sample_names.index(sample_name)
        return [value[index]
                for key, value in list(self.__key_to_sample_values.items())]

    def has_genotype_key(self, key):
        return key in self.__key_to_sample_values

    def genotype_keys(self):
        return list(self.__key_to_sample_values.keys())

    def get_sample_names(self):
        return self.__sample_names

    def has_sample(self, sample_name):
        return sample_name in self.__sample_names

    def to_vcf_header_columns(self):
        return ["FORMAT"] + self.__sample_names

    def to_vcf_columns(self):
        keys_string = [":".join(self.genotype_keys())]

        sample_strings = [[] for sample_name in self.__sample_names]

        for values_per_sample in list(self.__key_to_sample_values.values()):
            for index, value in enumerate(values_per_sample):
                if not isinstance(value, GenotypeCall) and value == []:
                    sample_strings[index].append('.')
                else:
                    sample_strings[index].append(to_vcf_str(value))

        return keys_string + [':'.join(sample_string)
                              for sample_string in sample_strings]

    def genotypes(self):
        genotypes = OrderedDict()
        for sample_name in self.__sample_names:
            genotypes[sample_name] = self.get_field(sample_name, GENOTYPE_KEY)
        return genotypes

    @property
    def has_merged_genotypes(self):
        return self.__merged_genotypes

    def merge_genotype_calls(self, dictionary_sample_name_to_genotype_call):
        # assert(not self.__merged_genotypes)
        for sample_name in self.__sample_names:
            sample_name_index = self.__sample_names.index(sample_name)
            self.__key_to_sample_values[GENOTYPE_KEY][sample_name_index] = merge_genotype_calls(
                self.__key_to_sample_values[GENOTYPE_KEY][sample_name_index],
                dictionary_sample_name_to_genotype_call[sample_name])
        self.__merged_genotypes = True


class SampleDataView(object):

    def __init__(self, sample_data, sample_name):
        self.__sample_data = sample_data
        self.__sample_name = sample_name

    def __contains__(self, key):
        return self.__sample_data.has_genotype_key(key)

    def __getitem__(self, key):
        return self.__sample_data.get_field(self.__sample_name, key)

    def keys(self):
        return self.__sample_data.genotype_keys()

    def values(self):
        return self.__sample_data.get_fields(self.__sample_name)

    def genotype(self):
        return self.__getitem__(GENOTYPE_KEY)
