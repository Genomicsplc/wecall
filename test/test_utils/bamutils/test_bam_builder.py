# All content Copyright (C) 2018 Genomics plc
import os

import pysam
from wecall.bamutils.bam_builder import BAMBuilder, RG_ID
from wecall.bamutils.sequence_bank import SequenceBank
from wecall.genomics.reference_chromosome import ReferenceChromosome
from wecall_test_drivers.base_test import BaseTest


class TestBamBuilder(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.chrom = "2"
        self.chrom_length = 200
        self.sample_name = "TEST_SAMPLE"
        self.filestub = "_"

    def test_header_for_multisample_multicontig(self):
        ref = ReferenceChromosome("")
        sequence_bank = SequenceBank(ref)
        builder = BAMBuilder(
            os.path.join(
                self.work_dir,
                self.filestub +
                ".bam"))
        builder.with_bam_contig_data("1", 10, "SAMPLE_ONE", sequence_bank)
        builder.with_bam_contig_data("2", 20, "SAMPLE_TWO", sequence_bank)

        expected_header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': 10, 'SN': "1"}, {'LN': 20, 'SN': "2"}], 'RG': [
            {"ID": RG_ID + "_SAMPLE_ONE", "SM": "SAMPLE_ONE"}, {"ID": RG_ID + "_SAMPLE_TWO", "SM": "SAMPLE_TWO"}]}

        self.assertDictEqual(expected_header, builder.header)

    def test_can_build_with_one_seq(self):
        ref = ReferenceChromosome("TCATAAAAAAAT")
        sequence_bank = SequenceBank(ref)
        sequence_bank.add_sequence(
            ".*G.........",
            "            ",
            n_fwd=2, n_rev=1
        )

        builder = BAMBuilder(
            os.path.join(
                self.work_dir,
                self.filestub +
                ".bam")) .with_bam_contig_data(
            self.chrom,
            self.chrom_length,
            self.sample_name,
            sequence_bank)
        builder.build()

        bam_file = pysam.Samfile(builder.filename, "rb")
        reads = list(bam_file.fetch())
        self.assertEqual(len(reads), 3)

        for read in reads:
            self.assertEqual(read.pos, 0)
            self.assertEqual(read.seq, "TGTAAAAAAAT")
            self.assertEqual(read.cigarstring, "1M1D10M")

        self.assertTrue(os.path.isfile(bam_file.filename))
        self.assertTrue(os.path.isfile(bam_file.filename.decode() + ".bai"))

    def test_can_build_with_defined_quality(self):
        ref = ReferenceChromosome("TCATAAAT")
        sequence_bank = SequenceBank(ref)
        sequence_bank.add_sequence(
            ".*G.....",
            "9 87  00",
            n_fwd=1, n_rev=0
        )

        builder = BAMBuilder(
            os.path.join(
                self.work_dir,
                self.filestub +
                ".bam")) .with_bam_contig_data(
            self.chrom,
            self.chrom_length,
            self.sample_name,
            sequence_bank)
        builder.build()

        bam_file = pysam.Samfile(builder.filename, "rb")
        reads = list(bam_file.fetch())
        self.assertEqual(len(reads), 1)
        self.assertEqual(reads[0].seq, "TGTAAAT")

        # ascii: "0": "!", "1": "+", "2": "5", "3": "?", "4": "H", "5": "S",
        # "6": "]", "7": "g", "8": "q", "9": "{"
        expected_qual = "{qgHH!!"
        self.assertEqual(reads[0].qual, expected_qual)

    def test_can_build_two_chroms(self):
        ref1 = ReferenceChromosome("TCATAAAAAAAT")
        sequence_bank1 = SequenceBank(ref1)
        sequence_bank1.add_sequence(".*G.........")

        ref2 = ReferenceChromosome("GGGG")
        sequence_bank2 = SequenceBank(ref2)
        sequence_bank2.add_sequence("..*.")

        builder = BAMBuilder(
            os.path.join(
                self.work_dir,
                self.filestub +
                ".bam")) .with_bam_contig_data(
            "1",
            100,
            "SAMPLE",
            sequence_bank1) .with_bam_contig_data(
                "X",
                50,
                "SAMPLE",
            sequence_bank2)
        builder.build()

        bam_file = pysam.Samfile(builder.filename, "rb")
        reads_chrom1 = list(bam_file.fetch(region="1:1-20"))
        self.assertEqual(len(reads_chrom1), 1)
        self.assertEqual(reads_chrom1[0].seq, "TGTAAAAAAAT")

        bam_file = pysam.Samfile(builder.filename, "rb")
        reads_chrom2 = list(bam_file.fetch(region="X:1-5"))
        self.assertEqual(len(reads_chrom2), 1)
        self.assertEqual(reads_chrom2[0].seq, "GGG")

        reads = list(bam_file.fetch())
        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].seq, "TGTAAAAAAAT")
        self.assertEqual(reads[1].seq, "GGG")

        self.assertRaises(ValueError, bam_file.fetch, region="2:1-20")

        self.assertTrue(os.path.isfile(bam_file.filename))
        self.assertTrue(os.path.isfile(bam_file.filename.decode() + ".bai"))
