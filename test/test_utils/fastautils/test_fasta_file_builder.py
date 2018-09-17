# All content Copyright (C) 2018 Genomics plc
import os
from unittest import expectedFailure

from wecall.fastautils.fasta_file_builder import FastaFileBuilder
from wecall_test_drivers.base_test import BaseTest
from pysam import Fastafile


class TestFastaFileBuilder(BaseTest):

    def __build_fasta_file(self, data):
        fasta_file = FastaFileBuilder(
            os.path.join(
                self.work_dir,
                "baa.fa"),
            line_length=5)
        for chrom, seq in list(data.items()):
            fasta_file.with_chrom(chrom, seq)
        return fasta_file.build().index()

    def test_should_be_able_to_fetch_section_of_genome(self):
        fasta_file = self.__build_fasta_file({'chr20': "TAGCATTATTATTATTATTATTATTA", })

        fasta_file = Fastafile(fasta_file.filename)
        self.assertEqual(fasta_file.fetch('chr20', 10, 20).upper(), "ATTATTATTA")

    def test_should_be_able_to_list_all_chromosomes(self):
        fasta_file = self.__build_fasta_file({'chr5': "T", 'chrX': "T", 'chr20': "T", })

        fasta_file = Fastafile(fasta_file.filename)
        self.assertEqual(sorted(fasta_file.references), sorted(['chr5', 'chr20', 'chrX']))

    def test_should_get_correct_chrom_length(self):
        chrom = 'chr20'
        seq = "TAGCATTATTATTATTATTATTATTA"
        fasta_file = self.__build_fasta_file({chrom: seq, })

        fasta_file = Fastafile(fasta_file.filename)
        self.assertEqual(fasta_file.get_reference_length(chrom), len(seq))

    @expectedFailure
    # "Cannot create reference files with dodgy chromosome content"
    def test_should_return_capitalised_sequence_from_ref_file(self):
        fasta_file = self.__build_fasta_file({'chr20': "tagcattattattattattattatta", })

        fasta_file = Fastafile(fasta_file.filename)
        self.assertEqual(fasta_file.fetch('chr20', 10, 20).upper(), "ATTATTATTA")
