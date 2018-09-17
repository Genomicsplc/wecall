# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestCallingVariantWithLowPercentageSupport(BaseTest):
    def test_somatic_call_on_homozygous_ref_background(self):
        chrom = "1"
        somatic_sample = "andy"

        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', chrom=chrom
        ).with_read(
            '........................................', n_rev=85, n_fwd=0, chrom=chrom, sample_name=somatic_sample
        ).with_read(
            '..................T.....................', n_rev=15, n_fwd=0, chrom=chrom, sample_name=somatic_sample
        ).with_ploidy(3)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(1) \
            .has_record_for_variant(Variant(chrom, 18, "C", "T")) \
            .with_sample(somatic_sample)

    def test_somatic_call_on_homozygous_alt_background(self):
        chrom = "1"
        somatic_sample = "andy"

        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', chrom=chrom
        ).with_read(
            '..........*.............................', n_rev=85, n_fwd=0, chrom=chrom, sample_name=somatic_sample
        ).with_read(
            '..........*.......T.....................', n_rev=15, n_fwd=0, chrom=chrom, sample_name=somatic_sample
        ).with_ploidy(3)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect \
            .has_record_for_variant(Variant(chrom, 18, "C", "T")) \
            .with_sample(somatic_sample)

    def test_somatic_call_on_het_alt_background(self):
        chrom = "1"
        somatic_sample = "andy"

        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            'CGGCGGTCGAACGGAGCCCCAAGCGAAGCTCAAAACATGG', chrom=chrom
        ).with_read(
            '..........T.............................', n_rev=50, n_fwd=0, chrom=chrom, sample_name=somatic_sample
        ).with_read(
            '..........*.............................', n_rev=35, n_fwd=0, chrom=chrom, sample_name=somatic_sample
        ).with_read(
            '..........*.......T.....................', n_rev=15, n_fwd=0, chrom=chrom, sample_name=somatic_sample
        ).with_ploidy(3)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect \
            .has_record_for_variant(Variant(chrom, 18, "C", "T")) \
            .with_sample(somatic_sample)
