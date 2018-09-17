# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.sequence_quality import SequenceQuality
import pysam
from wecall.bamutils.sample_bank import SampleBank
from wecall_test_drivers.ascii_wecall_runner import DEFAULT_SAMPLE_NAME
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.variant_caller_builder import VariantCallerBuilderFromSampleBank
import os


class AsciiQualityRecalibrationTest(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.sample_name = DEFAULT_SAMPLE_NAME
        self.output_stem = os.path.join(self.work_dir, "tmp")
        self.output_sam = "{}_{}.sam".format(
            self.output_stem, DEFAULT_SAMPLE_NAME)

    def assert_matching_ascii_qualities(self, first, second):
        non_ascii_first = SequenceQuality.parse_ascii_to_quality_string(first)
        non_ascii_second = SequenceQuality.parse_ascii_to_quality_string(
            second)
        self.assertEqual(non_ascii_first, non_ascii_second)

    def assert_quality_recalibrated_in_output_bam(
            self, ref_string, input_bam_seqs, output_bam_seqs):
        input_sample_bank = SampleBank(ref_string)
        input_sample_bank.add_sample_with_seqs_and_quals(
            self.sample_name, input_bam_seqs)

        output_sample_bank = SampleBank(ref_string)
        output_sample_bank.add_sample_with_seqs_and_quals(
            self.sample_name, output_bam_seqs)

        vc_builder = VariantCallerBuilderFromSampleBank(
            input_sample_bank, self.work_dir)
        vc_builder.configuration["recalibrateBaseQs"] = "true"
        vc_builder.configuration["intermediateRecalibFileStem"] = self.output_stem
        vc_builder.build().run()

        self.assertTrue(os.path.exists(self.output_sam))

        sam_file = pysam.Samfile(self.output_sam, "r")
        reads = list(sam_file.fetch())
        self.assertEqual(len(reads), len(output_sample_bank[self.sample_name]))

        # Sort the sam as in sequence bank.
        # output_sample_bank.sort_sequence_banks()
        output_reads = sorted(
            output_sample_bank[self.sample_name].build_reads(0, {}),
            key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq)
        )
        reads.sort(key=lambda x: (x.pos, x.seq, x.qual, x.cigarstring, x.mapq))

        for read, expected_sequence in zip(reads, output_reads):
            self.assertEqual(read.pos, expected_sequence.pos)
            self.assertEqual(read.seq, expected_sequence.seq)
            self.assert_matching_ascii_qualities(
                read.qual, expected_sequence.qual)
            self.assertEqual(read.cigarstring, expected_sequence.cigarstring)
            self.assertEqual(read.mapq, expected_sequence.mapq)

        sam_file.close()
