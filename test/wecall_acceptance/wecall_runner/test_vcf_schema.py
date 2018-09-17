# All content Copyright (C) 2018 Genomics plc
# -*- coding:utf8 -*-
from wecall_test_drivers.variant_caller_builder import VariantCallerBuilderFromSampleBank

from wecall.bamutils.sample_bank import SampleBank
from wecall.vcfutils.parser import VCFReaderContextManager
from wecall_test_drivers.base_test import BaseTest
import datetime
from wecall_test_drivers.wecall_schema import wecall_schema
import os


class TestVCFSchema(BaseTest):
    def __run_small_variant_caller(self, refcalls, format):
        sample_bank = SampleBank("T")
        sample_bank.add_sample_name("TEST").add_sequence(".")

        variant_caller_builder = VariantCallerBuilderFromSampleBank(sample_bank, self.work_dir)
        variant_caller_builder.configuration = {}  # clear config.
        variant_caller = variant_caller_builder.build()
        variant_caller.add_additional_command('outputRefCalls', refcalls)
        variant_caller.add_additional_command('outputFormat', "VCF{}".format(format))
        variant_caller.run()

        with VCFReaderContextManager(variant_caller.output_vcf) as vcf_file:
            actual_schema = vcf_file.read_header()

        reference = os.path.splitext(os.path.basename(
            variant_caller_builder.wecall_input_data.reference_filename))[0]
        expected_schema = wecall_schema(
            file_date=datetime.datetime.today().strftime('%F'),
            reference=reference,
            contigs={sample_bank.reference.chrom: {"length": sample_bank.reference.length_minus_deletions()}},
            add_ref_calls=refcalls,
            format=format)

        return expected_schema, actual_schema

    def test_correct_disclaimer(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertEqual(expected_schema.file_metadata['disclaimer'], actual_schema.file_metadata['disclaimer'])

    def test_correct_filedate_with_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertEqual(expected_schema.file_metadata['fileDate'], actual_schema.file_metadata['fileDate'])

    def test_correct_reference_with_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertEqual(expected_schema.file_metadata['reference'], actual_schema.file_metadata['reference'])

    def test_info_data_with_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertEqual(sorted(expected_schema.iter_info_data()), sorted(actual_schema.iter_info_data()))

    def test_info_data_without_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(False, "4.2")
        self.assertEqual(sorted(expected_schema.iter_info_data()), sorted(actual_schema.iter_info_data()))

    def test_sample_data_with_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertEqual(sorted(expected_schema.iter_sample_data()), sorted(actual_schema.iter_sample_data()))

    def test_correct_source_with_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertEqual(expected_schema.file_metadata['source'], actual_schema.file_metadata['source'])

    def test_we_call_outputs_contigs_with_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertEqual(sorted(expected_schema.iter_contigs()), sorted(actual_schema.iter_contigs()))

    def test_filters_schema_with_refcalls(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.2")
        self.assertCountEqual(sorted(expected_schema.iter_filters()), sorted(actual_schema.iter_filters()))

    def test_correct_filedate_with_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.1")
        self.assertEqual(expected_schema.file_metadata['fileDate'], actual_schema.file_metadata['fileDate'])

    def test_correct_reference_with_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.1")
        self.assertEqual(expected_schema.file_metadata['reference'], actual_schema.file_metadata['reference'])

    def test_info_data_with_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.1")
        self.assertEqual(sorted(expected_schema.iter_info_data()), sorted(actual_schema.iter_info_data()))

    def test_info_data_without_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(False, "4.1")
        self.assertEqual(sorted(expected_schema.iter_info_data()), sorted(actual_schema.iter_info_data()))

    def test_sample_data_with_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.1")
        self.assertEqual(sorted(expected_schema.iter_sample_data()), sorted(actual_schema.iter_sample_data()))

    def test_correct_source_with_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.1")
        self.assertEqual(expected_schema.file_metadata['source'], actual_schema.file_metadata['source'])

    def test_we_call_outputs_contigs_with_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.1")
        self.assertEqual(sorted(expected_schema.iter_contigs()), sorted(actual_schema.iter_contigs()))

    def test_filters_schema_with_refcalls_format_4_1(self):
        expected_schema, actual_schema = self.__run_small_variant_caller(True, "4.1")
        self.assertEqual(sorted(expected_schema.iter_filters()), sorted(actual_schema.iter_filters()))
