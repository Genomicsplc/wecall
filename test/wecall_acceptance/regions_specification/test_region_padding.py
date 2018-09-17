# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestRegionPadding(BaseTest):
    def test_should_call_snp_with_minimal_covering_region_using_default_padding(self):
        svc = SVCDriver(self)
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................T...........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:40-41')

        expect = svc.call()

        expect.with_output_vcf().record_count(1).has_record_for_variant(Variant('1', 40, "C", "T"))

    def test_should_call_del_with_minimal_covering_region_using_default_padding(self):
        svc = SVCDriver(self)
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................*...........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:40-41')

        expect = svc.call()

        expect.with_output_vcf().record_count(1).has_record_for_variant(Variant('1', 39, "GC", "G"))

    def test_should_call_del_with_minimal_covering_region_using_default_padding_with_region_before(self):
        svc = SVCDriver(self)
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................*...........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:39-40')

        expect = svc.call()

        expect.with_output_vcf().record_count(1).has_record_for_variant(Variant('1', 39, "GC", "G"))

    def test_should_not_call_del_if_region_doesnt_overlap_deleted_part(self):
        svc = SVCDriver(self)
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................*...........  ", chrom='1', n_fwd=10, n_rev=10)
        svc.with_region_string('1:38-39,1:41-42')

        expect = svc.call()

        expect.with_output_vcf().record_count(0)

    def test_should_call_first_snp_if_region_padding_is_zero(self):
        svc = SVCDriver(self)
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ...............................G..........  ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................T...........  ", chrom='1', n_fwd=6, n_rev=6)
        svc.with_region_string('1:0-41')
        svc.with_region_padding(0)

        expect = svc.call()

        expect.with_output_vcf().record_count(1).has_record_for_variant(Variant('1', 40, "C", "T"))

    def test_should_not_call_first_snp_if_region_padding_is_one(self):
        svc = SVCDriver(self)
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ...............................G..........  ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................T...........  ", chrom='1', n_fwd=6, n_rev=6)
        svc.with_region_string('1:0-41')
        svc.with_region_padding(1)

        expect = svc.call()

        expect.with_output_vcf().record_count(0)

    def test_should_cope_with_region_padding_which_pads_to_negative_index_into_reference(self):
        svc = SVCDriver(self)
        svc.with_ref_sequence(
            "AAAGCGTACAACCGGGTTAGTCACAAACCCGTTACGTATGCCCCCCCCCCCATG", chrom='1'
        ).with_read(
            "......................................................", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ...............................G..........  ", chrom='1', n_fwd=10, n_rev=10
        ).with_read(
            "          ..............................T...........  ", chrom='1', n_fwd=6, n_rev=6)
        svc.with_region_string('1:20-41')
        svc.with_region_padding(20)

        svc.call(expected_success=True)
