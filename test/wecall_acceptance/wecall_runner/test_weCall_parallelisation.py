# All content Copyright (C) 2018 Genomics plc
import os
import itertools
from wecall.bamutils.bam_builder import BAMBuilder
from wecall.bamutils.sequence_bank import SequenceBank
from wecall.fastautils.fasta_file_builder import FastaFileBuilder
from wecall.genomics.variant import Variant
from wecall.vcfutils.genotype_call import GenotypeCall
from wecall.wecall_utils.wecall_config_builder import WecallConfigBuilder
from wecall.wecall_utils.wecall_input_data import WecallInputData
from wecall_test_drivers.ascii_wecall_runner import DEFAULT_SAMPLE_NAME
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.variant_caller_wrapper import VariantCallerWrapper


class TestWeCallParallelisation(BaseTest):

    def setUp(self):
        BaseTest.setUp(self)
        self.chrom1 = "1"
        self.chrom2 = "2"
        #                   0         1            2         3         4
        # 012345678901   23456789012345678901234567890123456789
        self.ref_string1 = "AACCTTGGACGT***TATTCTGTCAATGCATCCCATTGCCGCCGCTAATCGCT"
        self.seq_string1 = " ...**......CTG.......*.....T...........A..T........."

        #                   0         1         2         3
        #                   0123456789012345678901234567890123456789
        self.ref_string2 = "GGGAATCATACATACTGGATTACCATTGGACCAGATTAGT"
        self.seq_string2 = "........*.............T................ "
        self.seq_string3 = "........*........ ....T................."
        self.block_size = 60

        self.sample_name1 = DEFAULT_SAMPLE_NAME
        self.sample_name2 = DEFAULT_SAMPLE_NAME

        self.vc_work_dir = os.path.join(self.work_dir, "vc_work_dir")
        os.makedirs(self.vc_work_dir)

    def setParallelAndSerialVariantCallers(self, copies1, copies2):
        '''Prepare the variant caller data for the test to run'''
        filestem = "vc_input"

        ref_file_builder = FastaFileBuilder(os.path.join(self.work_dir, filestem + ".fa"))
        ref1 = ref_file_builder.with_chrom(self.chrom1, self.ref_string1 * copies1)
        ref2 = ref_file_builder.with_chrom(self.chrom2, self.ref_string2 * copies2)

        self.repeat_length1 = ref1.length_minus_deletions() / copies1
        self.repeat_length2 = ref2.length_minus_deletions() / copies2

        ref_file_builder.build()
        ref_file_builder.index()

        seq_bank1 = SequenceBank(ref1)
        seq_bank1.add_sequence(self.seq_string1 * copies1, n_fwd=10, n_rev=10)

        seq_bank2 = SequenceBank(ref2)
        seq_bank2.add_sequence(self.seq_string2 * copies2, n_fwd=10, n_rev=10)
        seq_bank2.add_sequence(self.seq_string3 * copies2, n_fwd=10, n_rev=10)

        bam_builder = BAMBuilder(os.path.join(self.work_dir, filestem + ".bam"))
        bam_builder.with_bam_contig_data(ref1.chrom, ref1.length_minus_deletions(), self.sample_name1, seq_bank1)
        bam_builder.with_bam_contig_data(ref2.chrom, ref2.length_minus_deletions(), self.sample_name2, seq_bank2)
        bam_builder.build()

        wecall_input_data = WecallInputData([bam_builder.filename], ref_file_builder.filename)
        wecall_config_builder = WecallConfigBuilder(wecall_input_data, os.path.join(self.work_dir, filestem))
        wecall_config_builder.with_configuration("maxBlockSize", self.block_size)
        wecall_config_builder.with_configuration("noSimilarReadsFilter", False)
        wecall_config_builder.with_configuration("maxClusterDist", 20)
        wecall_config = wecall_config_builder.build()

        parallel_output_file_stem = os.path.join(self.work_dir, filestem + "_parallel")
        serial_output_file_stem = os.path.join(self.work_dir, filestem + "_serial")

        self.vc_wrapper_parallel = VariantCallerWrapper(parallel_output_file_stem, wecall_config)

        self.vc_wrapper_serial = VariantCallerWrapper(serial_output_file_stem, wecall_config)

    def test_should_give_same_results_in_parallel_as_in_series(self):
        self.setParallelAndSerialVariantCallers(1, 5)
        self.vc_wrapper_parallel.add_additional_command("numberOfJobs", "2")
        self.vc_wrapper_parallel.add_additional_command("workDir", self.vc_work_dir)
        self.vc_wrapper_parallel.run()
        self.vc_wrapper_serial.run()

        with open(self.vc_wrapper_parallel.output_vcf, "r") as parallel_vcf:
            with open(self.vc_wrapper_serial.output_vcf, "r") as serial_vcf:
                zipped = itertools.zip_longest(parallel_vcf, serial_vcf, fillvalue="MISSING_LINE")
                for parallel_vcf_line, serial_vcf_line in zipped:
                    if not parallel_vcf_line.startswith("##options"):
                        self.assertEqual(parallel_vcf_line, serial_vcf_line)

    def test_should_find_correct_variants(self):
        n_copies1 = 1
        n_copies2 = 5
        self.setParallelAndSerialVariantCallers(n_copies1, n_copies2)
        self.vc_wrapper_parallel.add_additional_command("numberOfJobs", "2")
        self.vc_wrapper_parallel.add_additional_command("workDir", self.vc_work_dir)
        self.vc_wrapper_parallel.add_additional_command("allowMNPCalls", False)
        self.vc_wrapper_parallel.run()

        expected_vars = set()
        for i in range(0, n_copies1):
            expected_vars.update({
                Variant(self.chrom1, 3 + i * self.repeat_length1, "CTT", "C"),
                Variant(self.chrom1, 11 + i * self.repeat_length1, "T", "TCTG"),
                Variant(self.chrom1, 18 + i * self.repeat_length1, "GT", "G"),
                Variant(self.chrom1, 25 + i * self.repeat_length1, "C", "T"),
                Variant(self.chrom1, 37 + i * self.repeat_length1, "G", "A"),
                Variant(self.chrom1, 40 + i * self.repeat_length1, "G", "T"),
            })

        for i in range(0, n_copies2):
            expected_vars.update({
                Variant(self.chrom2, 7 + i * self.repeat_length2, "AT", "A"),
                Variant(self.chrom2, 22 + i * self.repeat_length2, "C", "T"),
            })

        actual_parallel_variants = self.vc_wrapper_parallel.get_variant_callset(self).get_variants()
        self.assertEqual(expected_vars, actual_parallel_variants)

    def test_should_give_correct_output_for_different_sample_names(self):
        self.sample_name1 = "SAMPLE_A"
        self.sample_name2 = "SAMPLE_B"

        n_copies1 = 1
        n_copies2 = 5
        self.setParallelAndSerialVariantCallers(n_copies1, n_copies2)
        self.vc_wrapper_parallel.add_additional_command("numberOfJobs", "2")
        self.vc_wrapper_parallel.add_additional_command("workDir", self.vc_work_dir)
        self.vc_wrapper_parallel.run()

        expected_var_A_1 = Variant(self.chrom1, 3, "CTT", "C")
        expected_var_B_1 = Variant(self.chrom2, 7, "AT", "A")

        parallel_variants_with_genotypes = self.vc_wrapper_parallel \
            .get_variant_callset(self) \
            .get_variants_with_genotypes()

        self.assertTrue(expected_var_A_1 in list(parallel_variants_with_genotypes.keys()))
        self.assertTrue(expected_var_B_1 in list(parallel_variants_with_genotypes.keys()))

        self.assertEqual(GenotypeCall("1/1"), parallel_variants_with_genotypes[expected_var_A_1][self.sample_name1])
        self.assertEqual(GenotypeCall("./."), parallel_variants_with_genotypes[expected_var_A_1][self.sample_name2])
        self.assertEqual(GenotypeCall("./."), parallel_variants_with_genotypes[expected_var_B_1][self.sample_name1])
        self.assertEqual(GenotypeCall("1/1"), parallel_variants_with_genotypes[expected_var_B_1][self.sample_name2])
