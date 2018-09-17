# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.base_test import BaseTest
from wecall.bamutils.sample_bank import SampleBank
from wecall.wecall_utils.wecall_input_data_builder import WecallInputDataBuilder
import pysam


class TestWecallInputDataBuilder(BaseTest):
    def test_can_build_correct_ref_and_bam_file(self):
        bank = SampleBank("ATCCT*ATAATAAATAAATAAT")
        sample_name = "TEST_SAMPLE"
        bank.add_sample_name(sample_name)
        bank[sample_name].add_sequence("....CT.........T......")

        builder = WecallInputDataBuilder(self.work_dir).with_sample_bank(bank)

        input_files = builder.build()

        bam_file = pysam.Samfile(input_files.bam_filenames[0], "rb")
        for read in bam_file.fetch():
            self.assertEqual(read.pos, 0)
            self.assertEqual(read.seq, "ATCCCTATAATAAATTAATAAT")
            self.assertEqual(read.cigarstring, "5M1I16M")

        fasta_file = pysam.Fastafile(input_files.reference_filename)
        self.assertEqual(fasta_file.get_reference_length(bank.reference.chrom), 21)
        self.assertEqual(fasta_file.fetch(bank.reference.chrom, 0, 21), "ATCCTATAATAAATAAATAAT")

    def test_can_build_multiple_bam_files(self):
        bank = SampleBank("ATCCT*ATAATAAATAAATAAT")
        sample_name1 = "TEST_SAMPLE1"
        bank.add_sample_name(sample_name1)
        bank[sample_name1].add_sequence("....CT.........T......")

        sample_name2 = "TEST_SAMPLE2"
        bank.add_sample_name(sample_name2)
        bank[sample_name2].add_sequence(".....*.G..........*...")

        builder = WecallInputDataBuilder(self.work_dir).with_sample_bank(bank)
        input_bams = builder.build().bam_filenames

        bam_file1 = pysam.Samfile(input_bams[0], "rb")
        for read in bam_file1.fetch():
            self.assertEqual(read.pos, 0)
            self.assertEqual(read.seq, "ATCCCTATAATAAATTAATAAT")
            self.assertEqual(read.cigarstring, "5M1I16M")

        bam_file2 = pysam.Samfile(input_bams[1], "rb")
        for read in bam_file2.fetch():
            self.assertEqual(read.pos, 0)
            self.assertEqual(read.seq, "ATCCTAGAATAAATAAAAAT")
            self.assertEqual(read.cigarstring, "17M1D3M")

    def test_builds_correct_ref_and_bam_file_at_custom_position(self):
        pos_from = 100
        bank = SampleBank("ATCCT*ATAATAAATAAATAAT", pos_from)
        sample_name = "TEST_SAMPLE"
        bank.add_sample_name(sample_name)
        bank[sample_name].add_sequence("....CT.........T......")

        builder = WecallInputDataBuilder(self.work_dir).with_sample_bank(bank)

        input_files = builder.build()

        bam_file = pysam.Samfile(input_files.bam_filenames[0], "rb")
        for read in bam_file.fetch():
            self.assertEqual(read.pos, pos_from)
            self.assertEqual(read.seq, "ATCCCTATAATAAATTAATAAT")
            self.assertEqual(read.cigarstring, "5M1I16M")

        fasta_file = pysam.FastaFile(input_files.reference_filename)
        self.assertEqual(
            fasta_file.get_reference_length(bank.reference.chrom), 121)
        self.assertEqual(fasta_file.fetch(bank.reference.chrom, pos_from, pos_from + 21), "ATCCTATAATAAATAAATAAT")

    def test_should_be_able_to_build_bam_and_ref_data_with_multiple_chromosomes(self):
        bank_1 = SampleBank("A" * 10, 0, chrom='10')
        bank_1.add_sample_name("sample").add_sequence("." * 10)

        bank_2 = SampleBank("T" * 9, 0, chrom='20')
        bank_2.add_sample_name("sample").add_sequence("." * 9)

        builder = WecallInputDataBuilder(
            self.work_dir).with_sample_bank(bank_1).with_sample_bank(bank_2)

        input_files = builder.build()
        bam_file = pysam.Samfile(input_files.bam_filenames[0], "rb")

        for read in bam_file.fetch(reference='20'):
            self.assertEqual(read.pos, 0)
            self.assertEqual(read.seq, "T" * 9)
            self.assertEqual(read.cigarstring, "9M")

        for read in bam_file.fetch(reference='10'):
            self.assertEqual(read.pos, 0)
            self.assertEqual(read.seq, "A" * 10)
            self.assertEqual(read.cigarstring, "10M")
            print((dir(read)))
