# All content Copyright (C) 2018 Genomics plc
from wecall.wecall_utils.log_utils import log_timing_parser
import json

import os
from wecall_test_drivers.timed_command import TimedCommand
from wecall_test_drivers.tool_runner import log_file
from wecall_test_drivers.variant_callset import VariantCallSet

CANDIDATE_VARIANTS_FILE_KEY = "candidateVariantsFile"


class VariantCallerWrapper(object):

    def __init__(
            self,
            output_vcf_path_stem,
            wecall_config
    ):
        self.output_vcf = output_vcf_path_stem + ".vcf"
        self.log_filename = output_vcf_path_stem + ".log"
        self.wecall_config = wecall_config
        self.__additional_commands = {}
        self.__timmed_command = None

    @property
    def stderr(self):
        return self.__timmed_command.stderr

    @property
    def stdout(self):
        return self.__timmed_command.stdout

    @property
    def return_code(self):
        return self.__timmed_command.return_code

    @property
    def config_filename(self):
        return self.wecall_config.filename

    def get_variant_callset(self, test_case):
        variant_callset = VariantCallSet(test_case)
        variant_callset.add_vcf_variants(self.output_vcf)
        return variant_callset

    def add_additional_command(self, key, value):
        self.__additional_commands[key] = value

    def dump_timing_json(self, filename):
        times = self.__timmed_command.times
        with open(self.log_filename) as log_file:
            timing_data = log_timing_parser(log_file)
        per_file_timings = {}
        for timing_data_item in timing_data:
            assert(timing_data_item.timing_type == "IO")
            timed_file = timing_data_item.metadata["file"]
            if timed_file not in per_file_timings:
                per_file_timings[timed_file] = 0
            assert(timing_data_item.length_units == "us")
            per_file_timings[timed_file] += timing_data_item.length

        per_file_timings_in_seconds = {
            key: value / 10.0 ** 6 for key,
            value in per_file_timings.items()}
        times["IO"] = per_file_timings_in_seconds
        with open(filename, "w") as json_fp:
            json.dump(times, json_fp, indent=4, sort_keys=True)
            json_fp.write("\n")

    @property
    def system_time(self):
        return self.__timmed_command.system_time

    @property
    def user_time(self):
        return self.__timmed_command.user_time

    def run(self):
        cmd = [
            os.path.join(os.environ["WECALL_BIN"], "weCall"),
            "--output={}".format(self.output_vcf),
            "--logFilename={}".format(self.log_filename),
            "--config={}".format(self.wecall_config.filename)
        ]
        for key, value in self.__additional_commands.items():
            cmd.append("--{}={}".format(key, value))

        log_file(self.wecall_config.filename)

        # try:
        #     # remove file before trying to run tool
        #     os.unlink(self.output_vcf)
        # except OSError:
        #     pass

        self.__timmed_command = TimedCommand().start(cmd)

        log_file(self.output_vcf)

        return self
