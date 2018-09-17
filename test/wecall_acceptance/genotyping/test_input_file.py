# All content Copyright (C) 2018 Genomics plc
from os.path import join
from unittest import expectedFailure

from wecall.genomics.variant import Variant
from wecall.vcfutils.info_data import InfoData
from wecall.vcfutils.schema import Schema
from wecall.vcfutils.vcf_builder import VCFBuilder
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestInputSpecification(BaseTest):

    def test_doesnt_give_a_flying_damn_about_spurious_filter_header(self):
        chrom = "22"
        variant = Variant(chrom, 11, "A", "C")

        schema = Schema()
        complex_filter_name = '.+-*\\/~@?!%^&><=\"\'(){}[]_|'
        schema.set_filter(complex_filter_name, 'unusual characters')

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"), schema=schema)
        gv_builder.with_record_from_variant(variant, filters={complex_filter_name})
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

        expect = driver.call(expected_success=True)
        expect .with_output_vcf()\
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample)\
            .has_genotype("1/1")

    def test_doesnt_give_a_flying_damn_about_spurious_filters(self):
        chrom = "22"
        variant = Variant(chrom, 11, "A", "C")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(
            variant, filters={"#$.:@$%$%^&**()7!"})
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

        expect = driver.call(expected_success=True)
        expect.with_output_vcf()\
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample)\
            .has_genotype("1/1")

    def test_should_handle_complex_variant_input(self):
        chrom = "22"

        variant = Variant(chrom, 10, "CAA", "CA")

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
        expect.with_log()\
            .input_variant_trimmed_warning(variant, Variant(chrom, 11, "A", ""))
        expect.with_output_vcf()\
            .record_count(1)

    @expectedFailure  # "Unskip test if parameter made public"
    def test_should_raise_if_output_ref_calls_is_switched_on(self):
        chrom = "22"
        variant = Variant(chrom, 10, "CAA", "CA")

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
        ).with_output_ref_calls(True)

        driver.call(False).genotyping_is_incompatible_with_outputting_reference_calls_error()

    def test_doesnt_give_a_flying_damn_about_spurious_info(self):
        chrom = "22"
        variant = Variant(chrom, 11, "A", "C")

        gv_builder = VCFBuilder(join(self.work_dir, "genotype.vcf"))
        gv_builder.with_record_from_variant(variant,
                                            info=InfoData(None, {"#f$@$e%$%^&k**()7!": ["#o$@$f%$%f^&**()7!"]}))
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

        expect = driver.call(expected_success=True)
        expect.with_output_vcf() \
            .has_record_for_variant(variant)\
            .with_sample(dodgy_sample)\
            .has_genotype("1/1")
