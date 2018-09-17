# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestSmallRegions(BaseTest):

    def test_should_only_call_overlapping_ref_call(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_output_ref_calls(True).with_max_ref_call_size(1)

        svc.with_region_string('1:39-40')

        expect = svc.call()

        expect.with_output_vcf().record_count(1)

    def test_should_ignore_SNP_not_overlapping_region(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            ".......................................               ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................T...........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:0-40')

        expect = svc.call()

        expect.with_output_vcf().record_count(0)

    def test_should_not_ignore_SNP_not_overlapping_region(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            ".......................................               ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................T...........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:0-41')

        expect = svc.call()

        expect.with_output_vcf().record_count(1).has_record_for_variant(Variant('1', 40, "C", "T"))

    def test_should_use_read_evidence_outside_region_to_not_call_snp(self):
        svc = SVCDriver(self)
        # If no region was specified then only the SNP with stronger evidence
        # would be outputted.
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ...............................G..........  ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................T...........  ", chrom='1', n_fwd=6, n_rev=6)
        svc.with_region_string('1:0-41')

        expect = svc.call()

        expect.with_output_vcf().record_count(0)

    def test_should_ignore_deletion_not_overlapping_region(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            ".......................................               ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................****........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:0-39')

        expect = svc.call()

        expect.with_output_vcf().record_count(0)

    def test_should_include_deletion_overlapping_region(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            ".......................................               ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................****........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:0-41')

        expect = svc.call()

        expect.with_output_vcf().record_count(1).has_record_for_variant(Variant('1', 39, "GCCCC", "G"))

    def test_should_include_deletion_overlapping_two_region(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            ".......................................               ", chrom='1', n_fwd=10, n_rev=10).with_read(
            "          ..............................****........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:0-41,1:43-48')

        expect = svc.call()

        expect.with_output_vcf().record_count(1).has_record_for_variant(Variant('1', 39, "GCCCC", "G"))

    def test_should_call_variant_in_complex_region_within_small_calling_region(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTAGTCGGTAGGAATAATG", chrom='1').with_read(
            ".......................................               ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................****........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:0-44')

        expect = svc.call()

        expect.with_output_vcf().record_count(1)

    def test_should_not_ignore_variant_overlapping_edge_of_small_region(self):
        svc = SVCDriver(self)

        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            ".......................................               ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................****........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:0-42')

        expect = svc.call()

        expect.with_output_vcf().record_count(1)
