# All content Copyright (C) 2018 Genomics plc
import os
from wecall.bamutils.bam_builder import BAMBuilder, BAMContigData
from wecall.bamutils.sequence_bank import SequenceBank
from wecall.wecall_utils.wecall_input_data import WecallInputData
from wecall.fastautils.fasta_file_builder import FastaFileBuilder


class WecallInputDataBuilder(object):

    def __init__(self, work_dir):
        self.__sample_banks = None
        self.__work_dir = work_dir
        self.__ref_filename = None
        self.__bam_filenames = None

    def with_ref_filename(self, filename):
        self.__ref_filename = filename
        return self

    def with_bam_filenames(self, filenames):
        self.__bam_filenames = filenames
        return self

    def with_sample_bank(self, sample_bank):
        if self.__sample_banks is None:
            self.__sample_banks = []
        self.__sample_banks.append(sample_bank)
        return self

    def __build_ref(self):
        if self.__ref_filename is None:
            self.__ref_filename = os.path.join(self.__work_dir, "bah.fa")

        fasta_file_builder = FastaFileBuilder(self.__ref_filename)

        for sample_bank in self.__sample_banks:
            fasta_file_builder.with_chrom(
                sample_bank.reference.chrom,
                sample_bank.reference.fasta_string()
            )
        fasta_file_builder.build().index()
        return fasta_file_builder.filename

    def __build_bams(self):
        if len(self.__sample_banks) == 1:
            sample_names = self.__sample_banks[0].sample_names
        else:
            sample_names = set().union(
                *[sample_bank.sample_names for sample_bank in self.__sample_banks])

        bam_files = []
        if self.__bam_filenames is None:
            self.__bam_filenames = [
                os.path.join(
                    self.__work_dir,
                    "wecall_input_" +
                    sample_name +
                    ".bam") for sample_name in sample_names]

        for sample_name, filename in zip(sample_names, self.__bam_filenames):
            bam_file_builder = BAMBuilder(filename)
            for sample_bank in self.__sample_banks:
                sequence_bank = sample_bank.get(sample_name, None)
                if sequence_bank is not None:
                    bam_file_builder.with_bam_contig_data(
                        sample_bank.reference.chrom,
                        sample_bank.reference.length_minus_deletions(),
                        sample_name,
                        sequence_bank
                    )
            bam_file_builder.build()
            bam_files.append(bam_file_builder.filename)
        return bam_files

    def build(self):
        fasta_filename = self.__build_ref()
        bam_filenames = self.__build_bams()
        return WecallInputData(bam_filenames, fasta_filename)
