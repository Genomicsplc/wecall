# All content Copyright (C) 2018 Genomics plc
from wecall.fastautils.fasta_file_builder import FastaFileBuilder
from wecall_test_drivers.base_test import BaseTest
import os
import subprocess
from tempfile import NamedTemporaryFile
from wecall_test_drivers.svc_driver import SVCDriver


class TestCmdLineOptions(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)

        self.we_call = os.path.join(os.environ["ECHIDNA_BIN"], "weCall")
        self.bam_filename = os.path.join(self.work_dir, "input.bam")
        self.bam_index_filename = self.bam_filename + ".bai"
        self.ref_filename = os.path.join(self.work_dir, "refFile.fa")
        self.ref_index_filename = self.ref_filename + ".fai"
        self.output_filename = os.path.join(self.work_dir, "output.vcf")
        self.we_call_work_dir = os.path.join(self.work_dir, "temp_dir")
        self.log_filename = os.path.join(self.work_dir, "_.log")
        self.chrom = "2"
        self.chrom_string = "A" * 2

    def __build_default_fasta_file(self):
        fasta_file_builder = FastaFileBuilder(os.path.join(self.work_dir, "haha.fa"))
        fasta_file_builder.filename = self.ref_filename
        fasta_file_builder.with_chrom(self.chrom, self.chrom_string)
        fasta_file_builder.build().index()

    @property
    def default_cmd(self):
        return [
            self.we_call,
            "--inputs", self.bam_filename,
            "--refFile", self.ref_filename,
            "--output", self.output_filename,
            "--verbosity", "0",
            "--logFilename", self.log_filename,
        ]

    def test_should_not_predict_cmdline_option_for_input(self):
        self.__build_default_fasta_file()

        p = subprocess.Popen(
            [os.path.join(os.environ["ECHIDNA_BIN"], "weCall"), "--inp", "input/path"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()

        self.assertEqual(p.returncode, 1)
        self.assertEqual(stderr.decode(), "FAILED - unrecognised option '--inp'\n")

    def test_fail_when_input_files_missing(self):
        self.__build_default_fasta_file()

        p = subprocess.Popen(
            self.default_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()

        self.assertEqual(p.returncode, 1)
        self.assertRegex(
            stderr.decode(), "FAILED - File {} does not exist\n".format(self.bam_filename))

    def test_should_fail_if_input_file_missing_index_file(self):
        with NamedTemporaryFile(prefix="_bam", suffix=".bam", dir=self.work_dir) as fake_bam:
            self.bam_filename = fake_bam.name
            self.__build_default_fasta_file()

            p = subprocess.Popen(
                self.default_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()

            self.assertEqual(p.returncode, 1)
            self.assertRegex(
                stderr.decode(), "FAILED - Index file {}.bai does not exist\n".format(self.bam_filename))

    def test_should_fail_if_input_doesnt_have_bam_extension(self):
        with NamedTemporaryFile(prefix="_bam", suffix=".not_bam", dir=self.work_dir) as fake_bam:
            self.bam_filename = fake_bam.name

            self.__build_default_fasta_file()

            p = subprocess.Popen(
                self.default_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()

            self.assertEqual(p.returncode, 1)
            self.assertRegex(
                stderr.decode(), "FAILED - File {} does not have .bam extension\n".format(self.bam_filename)
            )

    def test_should_fail_if_reference_file_is_missing(self):
        open(self.bam_filename, "w").close()
        open(self.bam_index_filename, "w").close()

        p = subprocess.Popen(
            self.default_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()

        self.assertEqual(p.returncode, 1)
        self.assertRegex(stderr.decode(), "FAILED - File {} does not exist\n".format(self.ref_filename))

    def test_should_fail_if_reference_index_file_is_missing(self):
        open(self.bam_filename, "w").close()
        open(self.bam_index_filename, "w").close()

        with NamedTemporaryFile(prefix="_fa", suffix=".fa", dir=self.work_dir) as fake_fa:
            self.ref_filename = fake_fa.name
            p = subprocess.Popen(
                self.default_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()

            self.assertEqual(p.returncode, 1)
            self.assertRegex(stderr.decode(), "FAILED - Index file {}.fai does not exist\n".format(fake_fa.name))

    def test_should_fail_if_reference_doesnt_have_fa_extension(self):
        open(self.bam_filename, "w").close()
        open(self.bam_index_filename, "w").close()

        with NamedTemporaryFile(prefix="_fa", suffix=".not_fa", dir=self.work_dir) as fake_fa:
            self.ref_filename = fake_fa.name
            p = subprocess.Popen(
                self.default_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()

            self.assertEqual(p.returncode, 1)
            self.assertRegex(stderr.decode(), "FAILED - File {} does not have .fa extension\n".format(fake_fa.name))


class TestCmdLineOptionsUsingDriver(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.output_filename = os.path.join(self.work_dir, "output.vcf")
        self.we_call_work_dir = os.path.join(self.work_dir, "work_der")

    def test_should_fail_if_work_dir_is_not_directory_when_run_in_parallel(self):
        with NamedTemporaryFile(prefix="not_a_directory", dir=self.work_dir) as fake_dir:
            svc_driver = SVCDriver(self)
            svc_driver \
                .with_ref_sequence("A") \
                .with_read(".") \
                .with_number_of_jobs(2) \
                .with_work_dir(fake_dir.name) \
                .with_output_vcf_filename(self.output_filename)

            svc_driver\
                .call(expected_success=False)\
                .work_dir_not_a_directory_error(fake_dir.name)

    def test_should_fail_if_output_already_exists_serial(self):
        open(self.output_filename, "w").close()

        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_overwrite(False) \
            .with_output_vcf_filename(self.output_filename)

        assert(os.path.exists(self.output_filename))
        svc_driver.call(
            expected_success=False).output_exists_error(
            self.output_filename)

    def test_should_not_fail_due_to_existing_output_if_overwrite_specified(self):
        open(self.output_filename, "w").close()
        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_overwrite(True) \
            .with_output_vcf_filename(self.output_filename)
        svc_driver.call(expected_success=True)

    def test_should_not_fail_due_to_existing_output_if_no_option_specified(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence("A").with_read(".")
        svc_driver.call(expected_success=True)

    def test_should_fail_if_output_already_exists_parallel(self):
        open(self.output_filename, "w").close()

        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_overwrite(False) \
            .with_number_of_jobs(2) \
            .with_work_dir(self.we_call_work_dir) \
            .with_output_vcf_filename(self.output_filename)

        assert(os.path.exists(self.output_filename))
        svc_driver\
            .call(expected_success=False)\
            .output_exists_error(self.output_filename)

    def test_should_not_fail_due_to_existing_output_if_overwrite_specified_in_parallel(self):
        open(self.output_filename, "w").close()

        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_overwrite(True) \
            .with_number_of_jobs(2) \
            .with_work_dir(self.we_call_work_dir) \
            .with_output_vcf_filename(self.output_filename)

        assert(os.path.exists(self.output_filename))
        svc_driver.call(expected_success=True)

    def test_should_not_fail_due_to_existing_output_if_no_option_specified_in_parallel(self):
        open(self.output_filename, "w").close()

        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_number_of_jobs(2) \
            .with_work_dir(self.we_call_work_dir) \
            .with_output_vcf_filename(self.output_filename)

        assert(os.path.exists(self.output_filename))
        svc_driver.call(expected_success=True)

    def test_should_fail_if_mem_limit_is_negative(self):
        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_mem_limit(-100) \
            .with_output_vcf_filename(self.output_filename)

        svc_driver.call(False).with_mem_limit_range_error()

    def test_should_fail_if_mem_limit_is_just_below_acceptable_range(self):
        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_mem_limit(49) \
            .with_output_vcf_filename(self.output_filename)

        svc_driver.call(False).with_mem_limit_range_error()

    def test_should_fail_if_mem_limit_is_just_above_acceptable_range(self):
        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_mem_limit(1024 * 1024 + 1) \
            .with_output_vcf_filename(self.output_filename)

        svc_driver.call(False).with_mem_limit_range_error()

    def test_should_fail_if_invalid_output_file_format_provided(self):
        svc_driver = SVCDriver(self)
        svc_driver \
            .with_ref_sequence("A") \
            .with_read(".") \
            .with_output_format('bah4.1')

        svc_driver.call(False).with_incorrect_output_format_error()
