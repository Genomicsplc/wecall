# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestComplexIndelOutputRepresentations(BaseTest):
    def test_converts_del_plus_ins_to_snps(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACGCCCCTGC**AAAAAAAAAATCGTCTGTG", chrom=chrom
        ).with_read(
            "........**AT...................", n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            "..........**...................", n_fwd=10, n_rev=10, sample_name=sample_name)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(2)

        vcf_expect.has_record(chrom, 8, "G", "A")
        vcf_expect.has_record(chrom, 9, "C", "T")

    def test_always_puts_snp_on_right_1(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACGCCCCTGCAAAAAAAAAATCGTCTGTG", chrom=chrom
        ).with_read(
            "........T*...................", n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            ".............................", n_fwd=10, n_rev=10, sample_name=sample_name)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(2)

        vcf_expect.has_record(chrom, 7, "TG", "T")
        vcf_expect.has_record(chrom, 9, "C", "T")

    def test_always_puts_snp_on_right_2(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACGCCCCTGCAAAAAAAAAATCGTCTGTG", chrom=chrom
        ).with_read(
            "........*T...................", n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            ".............................", n_fwd=10, n_rev=10, sample_name=sample_name)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(2)

        vcf_expect.has_record(chrom, 7, "TG", "T")
        vcf_expect.has_record(chrom, 9, "C", "T")

    def test_joins_left_neighbouring_snp_to_deletion_and_simplifies_representation(self):
        sample_name = "sir_freedom_fries"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACGCCCCTGTAAAAAAAAAATCGTCTGTG", chrom=chrom
        ).with_read(
            "........T*...................", n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            ".............................", n_fwd=10, n_rev=10, sample_name=sample_name)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(1)
        vcf_expect.has_record(chrom, 7, "TG", "T")

    def test_switches_left_neighbouring_snp_to_insertion(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACGCCCCTG*AAAAAAAAAATCGTCTGTG", chrom=chrom
        ).with_read(
            "........TC...................", n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            ".........*...................", n_fwd=10, n_rev=10, sample_name=sample_name)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(2)

        vcf_expect.has_record(chrom, 6, "C", "CT")
        vcf_expect.has_record(chrom, 8, "G", "C")

    def test_leaves_insertion_on_left(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self).with_ref_sequence(
            "ACGCCCCT*CAAAAAAAAAATCGTCTGTG", chrom=chrom
        ).with_read(
            "........GT...................", n_fwd=10, n_rev=10, sample_name=sample_name
        ).with_read(
            "........*....................", n_fwd=10, n_rev=10, sample_name=sample_name)

        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.record_count(2)

        vcf_expect.has_record(chrom, 7, "T", "TG")
        vcf_expect.has_record(chrom, 8, "C", "T")
