# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from collections import OrderedDict


class WecallConfigFileRunnerTest(AsciiWecallRunnerTest):
    """
    Run a set of tests specified in a configuration file
    """

    def parse_file_into_list(self, input_file):
        config_data = []
        current_key = None

        for line in input_file:

            if line.strip().startswith("#"):
                continue

            line_no_newline = line.rstrip("\n")
            line_no_whitespace = line.strip()

            if not line_no_whitespace:
                continue

            if line_no_whitespace.endswith(":"):
                current_key = line_no_whitespace[0:-1]
                continue

            elif ":" in line:
                key, value = [x.strip() for x in line_no_whitespace.split(":")]
                config_data.append([key, [value]])
                current_key = None
                continue

            else:
                if len(config_data) == 0 or current_key is not None:
                    if current_key != config_data[-1][0]:
                        config_data.append([current_key, [line_no_newline]])
                    else:
                        config_data[-1][1].append(line_no_newline)
                else:
                    raise Exception("Invalid config file")

        return config_data

    def split_config_data(self, data):
        tests = OrderedDict()
        current_key = None

        for key, values in data:
            if key == "TestName":
                current_key = values[0]
                test_type = current_key.split("_")[0]
                tests[values[0]] = OrderedDict()
                tests[values[0]]["TestType"] = test_type
            else:
                if current_key is not None:
                    tests[current_key][key] = values
                else:
                    raise Exception("Invalid config file")
        return tests

    def check_variant_calling(self, test_name, test_data):
        samples = test_data["Samples"][0].split()
        reference_sequences = set()
        sample_bam_data = {}
        sample_variant_calls = {}

        for sample in samples:
            this_sample_bam_data = test_data["{}_Sequence".format(sample)]
            reference_sequence = this_sample_bam_data[0]
            reference_sequences.add(reference_sequence)
            reads = this_sample_bam_data[1:]
            sample_bam_data[sample] = reads
            expected_variants = test_data["{}_ExpectedVariants".format(sample)]

            for variant in expected_variants:
                chrom, pos, ref, alt, genotype = variant.split()
                this_variant = (int(pos), ref, alt)

                if this_variant not in sample_variant_calls:
                    sample_variant_calls[this_variant] = {}

                sample_variant_calls[this_variant][sample] = genotype

            assert len(reference_sequences) == 1
            the_reference = reference_sequences.pop()

        self.calls_variants_from_samples(
            the_reference,
            sample_bam_data,
            expected_call_stubs=sample_variant_calls
        )

    def check_recalibration(self, test_name, test_data):
        pass

    def run_from_config_file(self, fileName):
        """
        Load config file and run all specified tests
        """
        with open(fileName, 'r') as input_file:
            data = self.parse_file_into_list(input_file)
            test_data = self.split_config_data(data)

            for test in test_data:
                this_test_data = test_data[test]
                test_type = test.split("_")[0]

                if test_type == "VariantCalling":
                    self.check_variant_calling(test, this_test_data)
                elif test_type == "Recalibration":
                    self.check_recalibration(test, this_test_data)
                else:
                    raise Exception("Unknown test type {}".format(test_type))
