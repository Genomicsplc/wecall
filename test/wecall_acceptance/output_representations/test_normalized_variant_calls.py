# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestNormalizeVariantCallsOption(BaseTest):
    def test_should_normalize_to_having_large_indel_only_a_few_other_variants(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "GCCCCAGCCTCCCAAAGTGCATTGATTTTGTTGTTGTTGTGCTTATTTGCACTCCAGCCTGGCCTCTCCTTTCTTG", chrom=chrom
        ).with_read(
            "....................***********TGGGATTACAAGTGTGAACCATCGT....................",
            n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            "....................*****************.......................................",
            n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_normalize_variant_calls(True)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(9)

        vcf_expect.has_record(chrom, 19, "CATTGATTTTGTTGTTGT", "C").with_sample(
            sample_name).has_genotype("1|1")
        vcf_expect.has_record(chrom, 39, "T", "G").with_sample(
            sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 41, "C", "A").with_sample(
            sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 44, "A", "ACAAG").with_sample(
            sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 46, "T", "G").with_sample(
            sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 49, "C", "A").with_sample(
            sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 51, "C", "CCA").with_sample(
            sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 54, "C", "G").with_sample(
            sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 55, "A", "T").with_sample(
            sample_name).has_genotype("0|1")

    @expectedFailure
    def test_should_call_long_hom_deletion_separately_from_short_het_deletions(self):
        sample_name = "NA12878"
        chrom = "3"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACTACTTGTCCTTTCCGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCA"
            "TCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAAAATACAAAAAATTAGCCGGGCGTGGTAGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGA"
            "GGCAGGGGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAA"
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAACTACTTGTCCTTTCCTTTGTGCGTGTGTGCGTGTGTGTGTGTGTGT",
            chrom=chrom
        ).with_read(
            "ACTACTTGTCCTTTCC*************************************************************************************"
            "*****************************************************************************************************"
            "*****************************************************************************************************"
            "**************************************************************TTTGTGCGT******GTGTGTGTGTGTGTGT",
            n_fwd=20, n_rev=20, sample_name=sample_name, chrom=chrom
        ).with_read(
            "ACTACTTGTCCTTTCC*************************************************************************************"
            "*****************************************************************************************************"
            "*****************************************************************************************************"
            "**************************************************************TTTGTGCGTGT****GTGTGTGTGTGTGTGT",
            n_fwd=20, n_rev=20, sample_name=sample_name, chrom=chrom
        ).with_normalize_variant_calls(True)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(3)

        vcf_expect.has_record(
            chrom, 16,
            "CGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAG"
            "ATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAAAATACAAAAAATTAGCCGGGCGTGGTAG"
            "CGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGGGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGT"
            "GAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAA"
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAGAACTACTTGTCCTTTCC",
            "C"
        ).with_sample(sample_name).has_genotype("1|1")

        vcf_expect.has_record(chrom, 374, "TGTGTGC", "T").with_sample(sample_name).has_genotype("0|1")
        vcf_expect.has_record(chrom, 376, "TGTGC", "T").with_sample(sample_name).has_genotype("1|0")

    def test_should_call_basic_snps(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "GCCCCAGCCTCCCAAAGTGCATTGATTTTGTTGTTGTTGTGCTTATTTGCACTCCAGCCTGGCCTCTCCTTTCTTG", chrom=chrom
        ).with_read(
            "...............T.........A...............G..................................",
            n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            ".........................A..........................A.......................",
            n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_normalize_variant_calls(True)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(4)

        vcf_expect.has_record(chrom, 15, "A", "T").with_sample(
            sample_name).has_genotype("1|0")
        vcf_expect.has_record(chrom, 25, "T", "A").with_sample(
            sample_name).has_genotype("1|1")
        vcf_expect.has_record(chrom, 41, "C", "G").with_sample(
            sample_name).has_genotype("1|0")
        vcf_expect.has_record(chrom, 52, "T", "A").with_sample(
            sample_name).has_genotype("0|1")
