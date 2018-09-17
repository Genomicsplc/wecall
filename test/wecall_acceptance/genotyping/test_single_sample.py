# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall.genomics.variant import Variant
from wecall.vcfutils.vcf_builder import VCFBuilder
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver
from os.path import join


class TestSingleSampleGenotyping(BaseTest):
    def test_genotypes_variant_correctly_with_zero_coverage(self):
        chrom = "22"
        variant = Variant(chrom, 11, "A", "C")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().index()

        driver = SVCDriver(self)

        dodgy_sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA", chrom=chrom, pos_from=0
        ).with_read(
            "...........C.........", n_fwd=0, n_rev=0, chrom=chrom, sample_name=dodgy_sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename
        )

        expect = driver.call()
        expect.with_output_vcf()\
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample)\
            .has_genotype("./.")

    def test_genotypes_variant_correctly_complex_indel_which_is_snp_and_deletion(self):
        chrom = "22"
        variant = Variant(chrom, 10, "CA", "T")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().index()

        driver = SVCDriver(self)

        sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAAA", chrom=chrom, pos_from=0
        ).with_read(
            "..........T*..........", n_fwd=5, n_rev=5, chrom=chrom, sample_name=sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename
        )

        expect = driver.call()
        expect.with_output_vcf()\
            .has_record_for_variant(variant)\
            .with_sample(sample)\
            .has_genotype("1/1")

    def test_genotypes_variant_correctly_with_supporting_reads(self):
        chrom = "22"
        variant = Variant(chrom, 11, "A", "C")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().index()

        driver = SVCDriver(self)

        dodgy_sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA", chrom=chrom, pos_from=0
        ).with_read(
            "...........C.........", n_fwd=5, n_rev=5, chrom=chrom, sample_name=dodgy_sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename
        )

        expect = driver.call()
        expect.with_output_vcf()\
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample)\
            .has_genotype("1/1")

    def test_genotypes_mnp_correctly_with_supporting_reads(self):
        chrom = "22"
        variant = Variant(chrom, 11, "AAA", "CAC")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().index()

        driver = SVCDriver(self)

        dodgy_sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA", chrom=chrom, pos_from=0
        ).with_read(
            "...........C.C.......", n_fwd=5, n_rev=5, chrom=chrom, sample_name=dodgy_sample
        ).with_read(
            "...........C.........", n_fwd=5, n_rev=5, chrom=chrom, sample_name=dodgy_sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename
        )

        expect = driver.call()
        expect.with_output_vcf()\
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample)\
            .has_genotype("./1")

    @expectedFailure
    def test_gets_correct_genotype_if_not_fully_left_aligned(self):
        chrom = "22"

        variant = Variant(chrom, 12, "AA", "A")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().index()

        driver = SVCDriver(self)

        dodgy_sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAATACGCCCCCTACGCCCCCT", chrom=chrom, pos_from=0
        ).with_read(
            "...................*...................", n_fwd=10, n_rev=10, chrom=chrom, sample_name=dodgy_sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename
        )

        expect = driver.call()
        expect.with_output_vcf()\
            .record_count(1)\
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample).has_genotype("1/1")

    def test_genotypes_variant_correctly_with_unsupportive_reads(self):
        chrom = "22"
        variant = Variant(chrom, 11, "A", "C")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().index()

        driver = SVCDriver(self)

        dodgy_sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA", chrom=chrom, pos_from=0
        ).with_read(
            ".....................", n_fwd=5, n_rev=5, chrom=chrom, sample_name=dodgy_sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename
        )

        expect = driver.call()
        expect\
            .with_output_vcf()\
            .record_count(1) \
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample)\
            .has_genotype("0/0")

    def test_doesnt_call_extra_variants(self):
        chrom = "22"
        variant = Variant("1", 11, "A", "C")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().index()

        driver = SVCDriver(self)

        dodgy_sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA", chrom=chrom, pos_from=0
        ).with_read(
            "....G................", n_fwd=5, n_rev=5, chrom=chrom, sample_name=dodgy_sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename)

        expect = driver.call()
        expect \
            .with_output_vcf() \
            .record_count(0)

    def test_raises_if_genotyping_file_doesnt_exist(self):
        missing_file = join(
            self.work_dir,
            "I_DONT_EXIST_NOT_JUST_IN_THE_PHILOSOPHICAL_SENSE.vcf.gz")

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA",
        ).with_read(
            "....G................", n_fwd=5, n_rev=5,
        ).with_genotype_alleles(
            missing_file
        ).with_verbosity(0)

        driver.call(expected_success=False)\
            .missing_genotype_file(missing_file)

    def test_raises_if_genotyping_file_is_not_expected_format(self):

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.build()  # note: not compressed or indexed

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA",
        ).with_read(
            "....G................", n_fwd=5, n_rev=5,
        ).with_genotype_alleles(
            gv_builder.filename
        ).with_verbosity(0)

        driver.call(expected_success=False)\
            .unexpected_genotype_file_format(gv_builder.filename)

    def test_raises_if_genotyping_file_not_indexed(self):
        chrom = "22"
        variant = Variant(chrom, 11, "A", "C")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant)
        gv_builder.build().bgzip()  # note: no index

        driver = SVCDriver(self)

        dodgy_sample = "bobs_your_uncle"
        driver.with_ref_sequence(
            "ACGCCCCCTGCAAAAAAAAAA", chrom=chrom, pos_from=0
        ).with_read(
            "...........C.........", n_fwd=5, n_rev=5, chrom=chrom, sample_name=dodgy_sample
        ).with_genotype_alleles(
            gv_builder.compressed_filename
        ).with_verbosity(0)

        driver.call(expected_success=False)\
            .missing_genotype_index_file(gv_builder.compressed_filename_index)
