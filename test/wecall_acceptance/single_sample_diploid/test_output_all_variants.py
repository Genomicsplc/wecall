# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestOutputAllVariants(BaseTest):

    def test_should_output_isolated_variants_in_identical_way(self):
        ref_sequence = 'GTGG**AGACCTGAGCGAACAAGAGCGCAC'
        var_sequence = '  ..GA.........T.......**...  '

        normal = SVCDriver(self)
        normal \
            .with_ref_sequence(ref_sequence) \
            .with_read(var_sequence, n_fwd=10, n_rev=10)

        normal_vcf = normal \
            .call() \
            .with_output_vcf() \
            .record_count(3)

        all_variants = SVCDriver(self) \
            .with_all_variants(True)
        all_variants \
            .with_ref_sequence(ref_sequence) \
            .with_read(var_sequence, n_fwd=10, n_rev=10)

        all_variants_vcf = all_variants \
            .call() \
            .with_output_vcf() \
            .record_count(3)

        self.assertEqual(normal_vcf, all_variants_vcf)

    def test_should_output_each_individual_variant_from_MNP(self):
        chrom = 'chr12'
        driver = SVCDriver(self) \
            .with_all_variants(True) \
            .with_allow_MNP_calls(True)

        driver .with_ref_sequence(
            'GTGGGAAGACCTGAGCGAACAAGAGCGCAC', chrom=chrom
        ).with_read(
            '  .....T.....T.....T........  ', n_fwd=10, n_rev=10, chrom=chrom)

        vcf = driver.call().with_output_vcf()

        vcf.has_record_for_variant(Variant(chrom, 7, 'G', 'T')).with_filters({'LQ', 'NC'})
        vcf.has_record_for_variant(Variant(chrom, 13, 'A', 'T')).with_filters({'LQ', 'NC'})
        vcf.has_record_for_variant(Variant(chrom, 19, 'C', 'T')).with_filters({'LQ', 'NC'})
        vcf.has_record_for_variant(Variant(chrom, 7, 'GACCTGAGCGAAC', 'TACCTGTGCGAAT')).with_no_filters()
