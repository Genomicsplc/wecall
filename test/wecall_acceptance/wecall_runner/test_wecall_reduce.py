# All content Copyright (C) 2018 Genomics plc
from wecall.vcfutils.parser import VCFReaderContextManager
from wecall_test_drivers.base_test import BaseTest
from shutil import rmtree
from wecall.bamutils.sample_bank import SampleBank
from wecall_test_drivers.ascii_wecall_runner import DEFAULT_SAMPLE_NAME
from wecall_test_drivers.tool_runner import ToolRunner
from wecall_test_drivers.variant_caller_builder import VariantCallerBuilderFromSampleBank
from os import path, environ, mkdir
from subprocess import Popen, PIPE
import tempfile


class TestWeCallReduceCmdLine(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.tool_location = path.join(environ["ECHIDNA_BIN"], "weCall")

        self.input_directory_location = path.join(self.work_dir, "input_directory")
        self.log_filename = path.join(self.work_dir, "_.log")
        self.final_vcf_location = path.join(self.work_dir, "_.vcf")
        self.assertFalse(path.exists(self.final_vcf_location))

    def test_should_fail_if_input_directory_doesnt_exist(self):
        self.assertFalse(path.exists(self.input_directory_location))

        tool_runner = ToolRunner().start([
            self.tool_location,
            "reduce",
            "--inputDir={}".format(self.input_directory_location),
            "--output={}".format(self.final_vcf_location),
            "--verbosity={!s}".format(0),
            "--logFilename={!s}".format(self.log_filename),
        ])

        self.assertEqual(tool_runner.return_code, 1)
        self.assertRegex(tool_runner.stderr.decode(), "FAILED - Working dir: .* does not exist")

    def test_should_fail_if_directory_contains_files_with_non_vcf_extension(self):
        if path.exists(self.input_directory_location):
            rmtree(self.input_directory_location)
        mkdir(self.input_directory_location)

        self.assertTrue(path.exists(self.input_directory_location))
        tool_runner = ToolRunner()
        with tempfile.NamedTemporaryFile(prefix="fake_vcf", dir=self.input_directory_location):
            tool_runner.start([
                self.tool_location, "reduce",
                "--inputDir={}".format(self.input_directory_location),
                "--output={}".format(self.final_vcf_location),
                "--logFilename={!s}".format(self.log_filename),
            ])

        self.assertEqual(tool_runner.return_code, 1)
        self.assertRegex(tool_runner.stderr.decode(), "FAILED - file .* is not a VCF")

    def test_should_fail_if_directory_is_empty(self):
        if path.exists(self.input_directory_location):
            rmtree(self.input_directory_location)
        mkdir(self.input_directory_location)

        self.assertTrue(path.exists(self.input_directory_location))
        tool_runner = ToolRunner() \
            .start([
                self.tool_location, "reduce",
                "--inputDir={}".format(self.input_directory_location),
                "--output={}".format(self.final_vcf_location),
                "--logFilename={!s}".format(self.log_filename),
            ])

        self.assertEqual(tool_runner.return_code, 1)
        self.assertRegex(tool_runner.stderr.decode(),
                         "FAILED - directory .*{} is empty".format(self.input_directory_location))

    def test_should_fail_if_directory_contains_a_directory(self):
        if path.exists(self.input_directory_location):
            rmtree(self.input_directory_location)
        mkdir(self.input_directory_location)

        self.assertTrue(path.exists(self.input_directory_location))

        non_vcf_file = tempfile.mkdtemp(
            prefix="fake_file",
            suffix=".vcf",
            dir=self.input_directory_location)

        tool_runner = ToolRunner().start([
            self.tool_location, "reduce",
            "--inputDir={}".format(self.input_directory_location),
            "--output={}".format(self.final_vcf_location),
            "--logFilename={}".format(self.log_filename),
        ])

        self.assertEqual(tool_runner.return_code, 1)
        self.assertRegex(
            tool_runner.stderr.decode(),
            "FAILED - .*{} is not a file".format(non_vcf_file))
        if path.exists(non_vcf_file):
            rmtree(non_vcf_file)

    def test_should_fail_if_directory_contains_files_with_empty_vcfs(self):
        if path.exists(self.input_directory_location):
            rmtree(self.input_directory_location)
        mkdir(self.input_directory_location)

        self.assertTrue(path.exists(self.input_directory_location))
        with tempfile.NamedTemporaryFile(prefix="fake_vcf", suffix=".vcf", dir=self.input_directory_location):
            p = Popen(
                [
                    self.tool_location, "reduce",
                    "--inputDir={}".format(self.input_directory_location),
                    "--output={}".format(self.final_vcf_location),
                    "--logFilename={!s}".format(self.log_filename),
                ],
                stdout=PIPE,
                stderr=PIPE
            )
            stdout, stderr = p.communicate()

            self.assertEqual(p.returncode, 1)
            self.assertRegex(stderr.decode(), "FAILED - file .* is not a valid VCF")

    def test_should_fail_if_other_params_are_passed(self):
        if path.exists(self.input_directory_location):
            rmtree(self.input_directory_location)
        mkdir(self.input_directory_location)

        self.assertTrue(path.exists(self.input_directory_location))

        p = Popen(
            [
                self.tool_location, "reduce",
                "--inputDir={}".format(self.input_directory_location),
                "--output={}".format(self.final_vcf_location),
                "--regions=1",
                "--logFilename={!s}".format(self.log_filename),
            ],
            stdout=PIPE,
            stderr=PIPE
        )
        stdout, stderr = p.communicate()

        self.assertEqual(p.returncode, 1)
        self.assertEqual(
            stderr.decode(), "FAILED - unrecognised option '--regions=1'\n")


class TestReduce(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.tool_location = path.join(environ["ECHIDNA_BIN"], "weCall")
        self.intermediate_vcfs_dir = path.join(self.work_dir, "intermediate_vcfs")
        self.final_vcf_location = path.join(self.work_dir, "reduced.vcf")
        self.log_filename = path.join(self.work_dir, "_.log")

        mkdir(self.intermediate_vcfs_dir)

        self.stdout, self.stderr = None, None

    def tearDown(self):
        if path.exists(self.intermediate_vcfs_dir):
            rmtree(self.intermediate_vcfs_dir)
        BaseTest.tearDown(self)

    def test_should_reduce_a_wecall_produced_vcf_to_a_valid_vcf(self):
        temp_vcf_filename = self.__run_wecall_variant_caller(
            "1",
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["...................T......................"],
        )

        with VCFReaderContextManager(temp_vcf_filename) as temp_vcf:
            reference_header = temp_vcf.header
            reference_records = list(temp_vcf.read_records())

        self.__run_wecall_reduce()

        with VCFReaderContextManager(self.final_vcf_location) as final_vcf:
            self.assertEqual(final_vcf.header, reference_header)
            final_records = list(final_vcf.read_records())

            self.assertEqual(len(final_records), 1)
            self.assertEqual(final_records, reference_records)

    def test_should_derive_use_lexigraphical_order_of_vcfs_for_reduce(self):
        temp_vcf_filename_b = self.__run_wecall_variant_caller(
            "2",
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["............T........................C...."], vcf_stem="ab"
        )
        temp_vcf_filename_a = self.__run_wecall_variant_caller(
            "1",
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["...................T......................"], vcf_stem="aa"
        )

        with VCFReaderContextManager(temp_vcf_filename_a) as temp_vcf_a:
            with VCFReaderContextManager(temp_vcf_filename_b) as temp_vcf_b:
                # aa is lexicographical less than ab
                reference_records = list(
                    temp_vcf_a.read_records()) + list(temp_vcf_b.read_records())

        self.__run_wecall_reduce()

        with VCFReaderContextManager(self.final_vcf_location) as final_vcf:
            final_records = list(final_vcf.read_records())

            self.assertEqual(len(final_records), 3)
            self.assertEqual(final_records, reference_records)

    def test_should_obtain_correct_vcf_header_on_reduce(self):
        temp_vcf_filename_b = self.__run_wecall_variant_caller(
            "2",
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["............T........................C...."], vcf_stem="ab"
        )
        temp_vcf_filename_a = self.__run_wecall_variant_caller(
            "1",
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["...................T......................"], vcf_stem="aa"
        )

        with VCFReaderContextManager(temp_vcf_filename_a) as temp_vcf_a:
            with VCFReaderContextManager(temp_vcf_filename_b) as temp_vcf_b:
                temp_vcf_a.read_header()
                header_a = temp_vcf_a.header

                temp_vcf_b.read_header()
                header_b = temp_vcf_b.header

        self.__run_wecall_reduce()

        with VCFReaderContextManager(self.final_vcf_location) as final_vcf:
            final_vcf.read_header()

            expected_header = header_a

            expected_header.set_contig('2', header_b.get_contig('2').length)
            self.assertEqual(final_vcf.header, expected_header)

    def __run_wecall_variant_caller(self, chrom, reference_string, sequence_list, vcf_stem=None):
        if vcf_stem is None:
            vcf_stem = chrom
        sample_bank = SampleBank(reference_string, chrom=chrom)
        sample_bank.add_sample_with_seqs_and_quals(DEFAULT_SAMPLE_NAME, sequence_list, n_fwd=10, n_rev=10)
        vc_builder = VariantCallerBuilderFromSampleBank(sample_bank, self.work_dir)
        vc_wrapper = vc_builder.build()
        vc_wrapper.add_additional_command("allowMNPCalls", False)
        vc_wrapper.output_vcf = path.join(self.intermediate_vcfs_dir, "{}.vcf".format(vcf_stem))
        vc_wrapper.run()
        return vc_wrapper.output_vcf

    def __run_wecall_reduce(self):
        cmd = [
            self.tool_location,
            "reduce",
            "--inputDir={}".format(self.intermediate_vcfs_dir),
            "--output={}".format(self.final_vcf_location),
            "--logFilename={!s}".format(self.log_filename),
        ]

        print("Running '{}'".format(" ".join(cmd)))

        p = Popen(cmd, stdout=PIPE, stderr=PIPE)
        self.stdout, self.stderr = p.communicate()
        self.assertEqual(p.returncode, 0, "stdout = {}\nstderr={}".format(self.stdout, self.stderr))
        self.assertTrue(path.exists(self.final_vcf_location))
