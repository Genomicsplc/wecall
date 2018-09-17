# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.sample_bank import SampleBank
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall.vcfutils.genotype_call import GenotypeCall as GT
from wecall_test_drivers.variant_caller_builder import VariantCallerBuilderFromSampleBank
from wecall_test_drivers.vcf_expectation import VCFExpectation


class TestSNPGenotyping(AsciiWecallRunnerTest):

    def test_genotypes_snp_as_heterozygous_when_supported_by_2_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["..........................................",
             "..........................................",
             "..........................................",
             ".....................C....................",
             ".....................C...................."],

            [(21, "A", "C", {"GT": GT("1/0"),
                             "GQ": [66],
                             "DP": [5],
                             "AD": [3, 2]})]
        )

    def test_genotypes_snp_as_heterozygous_when_supported_by_3_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["..........................................",
             "..........................................",
             ".....................C....................",
             ".....................C....................",
             ".....................C...................."],

            [(21, "A", "C", {"GT": GT("1/0"),
                             "GQ": [66],
                             "DP": [5],
                             "AD": [2, 3]})]
        )

    def test_genotypes_snp_as_heterozygous_when_supported_by_4_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["..........................................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C...................."],

            [(21, "A", "C", {"GT": GT("1/0"),
                             "GQ": [27],
                             "DP": [5],
                             "AD": [1, 4]})]
        )

    def test_genotypes_snp_as_homozygous_when_supported_by_5_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            [".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C...................."],

            [(21, "A", "C", {"GT": GT("1/1"),
                             "GQ": [12],
                             "DP": [5],
                             "AD": [0, 5]})]
        )

    def test_genotypes_snp_as_heterozygous_when_supported_by_9_reads_out_of_10(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            [".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             "..........................................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C...................."],

            [(21, "A", "C", {"GT": GT("0/1"),
                             "GQ": [12],
                             "DP": [10],
                             "AD": [1, 9]})]
        )

    def test_genotypes_snp_as_homozygous_when_supported_by_15_reads_out_of_16(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            [".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             "..........................................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C....................",
             ".....................C...................."],

            [(21, "A", "C", {"GT": GT("1/1"),
                             "GQ": [7],
                             "DP": [16],
                             "AD": [1, 15]})]
        )


class TestDeletionGenotyping(AsciiWecallRunnerTest):

    def test_genotypes_one_base_deletion_as_heterozygous_when_supported_by_2_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["..........................................",
             "..........................................",
             "..........................................",
             ".....................*....................",
             ".....................*...................."],

            [(20, "CA", "C", {"GT": GT("1/0"),
                              "DP": [5],
                              "AD": [3, 2]})]
        )

    def test_genotypes_one_base_deletion_as_heterozygous_when_supported_by_3_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["..........................................",
             "..........................................",
             ".....................*....................",
             ".....................*....................",
             ".....................*...................."],

            [(20, "CA", "C", {"GT": GT("1/0"),
                              "DP": [5],
                              "AD": [2, 3]})]
        )

    def test_genotypes_one_base_deletion_as_heterozygous_when_supported_by_4_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["..........................................",
             ".....................*....................",
             ".....................*....................",
             ".....................*....................",
             ".....................*...................."],

            [(20, "CA", "C", {"GT": GT("1/0"),
                              "DP": [5],
                              "AD": [1, 4]})]
        )

    def test_genotypes_one_base_deletion_as_homozygous_when_supported_by_5_reads_out_of_5(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            [".....................*....................",
             ".....................*....................",
             ".....................*....................",
             ".....................*....................",
             ".....................*...................."],

            [(20, "CA", "C", {"GT": GT("1/1"),
                              "DP": [5],
                              "AD": [0, 5]})]
        )


class TestGenotypeLikelihoodValues(BaseTest):
    def test_should_have_near_zero_AA_genotype_likelihood_for_hom_alt_call(self):
        chr1 = 'chr1'
        sample_bank = SampleBank("TTTTTAAAAAAAAAAAAAAAAAAAA", chrom=chr1)

        sequence_bank = sample_bank.add_sample_name('sample_1')

        sequence_bank.add_sequence(
            "............C............", n_fwd=20, n_rev=20)

        vc_wrapper_builder = VariantCallerBuilderFromSampleBank(
            sample_bank, self.work_dir)
        variant_output = vc_wrapper_builder.build().run().output_vcf

        vcf_expectation = VCFExpectation(self, variant_output)
        record_expectation = vcf_expectation.has_record_for_variant(Variant(chr1, 12, "A", "C"))
        sample_expectation = record_expectation.with_sample("sample_1")

        sample_expectation.has_genotype("1|1").has_AA_genotype_likelihood(0)

    def test_should_have_near_zero_RA_genotype_likelihood_for_het_call(self):
        chr1 = 'chr1'
        sample_bank = SampleBank("TTTTTAAAAAAAAAAAAAAAAAAAA", chrom=chr1)

        sequence_bank = sample_bank.add_sample_name('sample_1')

        sequence_bank.add_sequence(
            "............C............", n_fwd=20, n_rev=20)
        sequence_bank.add_sequence(
            ".........................", n_fwd=20, n_rev=20)

        vc_wrapper_builder = VariantCallerBuilderFromSampleBank(
            sample_bank, self.work_dir)
        variant_output = vc_wrapper_builder.build().run().output_vcf

        vcf_expectation = VCFExpectation(self, variant_output)
        record_expectation = vcf_expectation.has_record_for_variant(Variant(chr1, 12, "A", "C"))
        sample_expectation = record_expectation.with_sample("sample_1")

        sample_expectation.has_genotype("1|0").has_RA_genotype_likelihood(0)

    def test_should_have_near_zero_RR_genotype_likelihood_for_hom_alt_call(self):
        chr1 = 'chr1'
        sample_bank = SampleBank("TTTTTAAAAAAAAAAAAAAAAAAAA", chrom=chr1)

        sequence_bank_1 = sample_bank.add_sample_name('sample_1')
        sequence_bank_1.add_sequence(
            ".........................", n_fwd=20, n_rev=20)

        sequence_bank_2 = sample_bank.add_sample_name('sample_2')
        sequence_bank_2.add_sequence(
            "............C............", n_fwd=20, n_rev=20)

        vc_wrapper_builder = VariantCallerBuilderFromSampleBank(sample_bank, self.work_dir)
        variant_output = vc_wrapper_builder.build().run().output_vcf

        vcf_expectation = VCFExpectation(self, variant_output)
        record_expectation = vcf_expectation.has_record_for_variant(Variant(chr1, 12, "A", "C"))
        sample_expectation = record_expectation.with_sample("sample_1")

        sample_expectation.has_genotype("0|0").has_RR_genotype_likelihood(0.0)
