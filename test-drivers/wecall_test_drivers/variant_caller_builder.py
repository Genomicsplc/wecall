# All content Copyright (C) 2018 Genomics plc
from wecall.wecall_utils.wecall_input_data_builder import WecallInputDataBuilder
from wecall_test_drivers.variant_caller_wrapper import VariantCallerWrapper
import os

from wecall_test_drivers import variant_caller_wrapper

# The standard human reference genome. Needed for variant calling
from wecall.wecall_utils.wecall_config_builder import WecallConfigBuilder

EDNA_FATAL = 0
EDNA_ERROR = 1
EDNA_WARNING = 2
EDNA_INFO = 3
EDNA_TIMING = 4
EDNA_DEBUG = 5
EDNA_SUPER_DEBUG = 6


class VariantCallerBuilder(object):
    def __init__(
            self,
            workdir,
            dataset,
            chrom_intervals=None,
    ):
        self.__output_file_path_stem = os.path.join(workdir, "weCall")

        self.__wecall_config_builder = WecallConfigBuilder(
            dataset, self.__output_file_path_stem)

        if chrom_intervals is not None:
            self.with_configuration("regions", ",".join(
                [str(intv) for intv in chrom_intervals]))

        self.with_configuration("verbosity", EDNA_FATAL)
        self.with_configuration("overwrite", True)

    def with_configuration(self, key, value):
        if value is not None:
            self.__wecall_config_builder.with_configuration(key, value)
        return self

    def with_configuration_dict(self, config_dict):
        for key, value in config_dict.items():
            self.with_configuration(key, value)
        return self

    def build(self):
        wecall_config = self.__wecall_config_builder.build()
        return variant_caller_wrapper.VariantCallerWrapper(
            self.__output_file_path_stem, wecall_config)


class VariantCallerBuilderFromSampleBank(object):
    def __init__(self, sample_bank, work_dir):
        self.sample_bank = sample_bank
        self.work_dir = work_dir
        self.filestem = os.path.join(self.work_dir, "_")
        self.configuration = {'noSimilarReadsFilter': False, 'minCallQual': 2}
        self.wecall_input_data = None

    def build(self):
        self.wecall_input_data = WecallInputDataBuilder(
            self.work_dir).with_sample_bank(
            self.sample_bank).build()
        wecall_config_builder = WecallConfigBuilder(
            self.wecall_input_data, self.filestem)

        for key, value in self.configuration.items():
            wecall_config_builder.with_configuration(key, value)
        wecall_config = wecall_config_builder.build()

        return VariantCallerWrapper(self.filestem, wecall_config)
