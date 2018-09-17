# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver

ref_alt = "<NON_REF>"


class TestCallingWithPloidy1(BaseTest):
    def test_should_call_variants(self):
        chrom = 'chr1'
        sample_name = 'sample'
        svc = SVCDriver(self) \
            .with_ploidy(1)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTC***AACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "................G.....ATG.......***.........", n_rev=10, n_fwd=10, chrom=chrom, sample_name=sample_name
        )

        expect = svc.call()

        vcf = expect \
            .with_output_vcf() \
            .record_count(3)

        vcf.has_record_for_variant(Variant(chrom, 16, 'T', 'G')).with_sample(sample_name).has_genotype('1')
        vcf.has_record_for_variant(Variant(chrom, 21, 'C', 'CATG')).with_sample(sample_name).has_genotype('1')
        vcf.has_record_for_variant(Variant(chrom, 28, 'TTAC', 'T')).with_sample(sample_name).has_genotype('1')

    def test_should_support_refcalls(self):
        chrom = 'chr1'
        sample_name = 'sample'
        svc = SVCDriver(self) \
            .with_ploidy(1) \
            .with_output_ref_calls(True)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCTCAAACCCGTTACGTATGCATG", chrom=chrom
        ).with_read(
            "............................................", n_rev=10, n_fwd=10, chrom=chrom, sample_name=sample_name
        )

        expect = svc.call()

        vcf = expect \
            .with_output_vcf() \
            .record_count(1)

        vcf.has_record_for_variant(Variant(chrom, 0, 'A', ref_alt)).with_sample(sample_name).has_genotype('0')
