# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestInfoAnnotations(BaseTest):
    def test_should_have_standard_set_of_info_keys_for_variant(self):
        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "TAGAATTGTTTGAGCTCTTTGTATTTCCTGTTATTAATCCCTTGTCAGAAGGGTCGTTTG", )\
            .with_read(
                "....................A.......................................", n_fwd=3, n_rev=0)

        svc_driver.call()\
            .with_output_vcf()\
            .record_count(1)\
            .has_record_for_variant(Variant("1", 20, "G", "A"))\
            .with_info()\
            .has_keys("PP", "DP", "DPR", "DPF", "VC", "VCR", "VCF", "ABPV", "SBPV", "MQ", "QD", "BR")

    def test_snps_at_same_location_should_have_expected_coverage(self):
        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "TAGAATTGTTTGAGCTCTTTGTATTTCCTGTTATTAATCCCTTGTCAGAAGGGTCGTTTG",)\
            .with_read(
                "....................A.......................................", n_fwd=27, n_rev=12)\
            .with_read(
                "....................T.......................................", n_fwd=0, n_rev=16)

        vcf_expect = svc_driver.call().with_output_vcf().record_count(2)

        vcf_expect\
            .has_record_for_variant(Variant("1", 20, "G", "A"))\
            .with_info().with_field("DP", [55]).with_field("VC", [39])

        vcf_expect\
            .has_record_for_variant(Variant("1", 20, "G", "T"))\
            .with_info().with_field("DP", [55]).with_field("VC", [16])

    def test_snp_and_insertion_at_same_location_should_have_expected_coverage(self):
        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "TAGAATTGTTTGAGCTCTTTG**TATTTCCTGTTATTAATCCCTTGTCAGAAGGGTCGTTTG",)\
            .with_read(
                "....................A**.......................................", n_fwd=27, n_rev=12)\
            .with_read(
                ".....................AT.......................................", n_fwd=0, n_rev=16)

        vcf_expect = svc_driver.call().with_output_vcf().record_count(2)

        vcf_expect\
            .has_record_for_variant(Variant("1", 20, "G", "A"))\
            .with_info().with_field("DP", [55]).with_field("VC", [39])

        vcf_expect\
            .has_record_for_variant(Variant("1", 20, "G", "GAT"))\
            .with_info().with_field("DP", [55]).with_field("VC", [16])

    def test_snp_and_deletion_at_same_location_should_have_expected_coverage(self):
        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "TAGAATTGTTTGAGCTCTTTGTATTTCCTGTTATTAATCCCTTGTCAGAAGGGTCGTTTG",)\
            .with_read(
                "....................A.......................................", n_fwd=27, n_rev=12)\
            .with_read(
                ".....................**.....................................", n_fwd=0, n_rev=16)

        vcf_expect = svc_driver.call().with_output_vcf().record_count(2)

        vcf_expect\
            .has_record_for_variant(Variant("1", 20, "G", "A"))\
            .with_info().with_field("DP", [55]).with_field("VC", [39])

        vcf_expect\
            .has_record_for_variant(Variant("1", 20, "GTA", "G"))\
            .with_info().with_field("DP", [55]).with_field("VC", [16])
