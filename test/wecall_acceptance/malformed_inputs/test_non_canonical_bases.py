# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestNonCanincalBases(BaseTest):
    def test_should_not_output_snp_with_unknown_character_in_ref(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATNGATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "..............T.........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)

    def test_should_not_output_snp_with_unknown_character_in_alt(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATGGATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "..............N.........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)

    def test_should_not_output_snp_with_unknown_base_before_insertion(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATN*ATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "...............T........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)

    def test_should_not_output_snp_with_unknown_base_before_deletion(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATNTATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "...............*........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)

    def test_should_not_output_snp_with_non_canonical_character_in_ref(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATKGATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "..............T.........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)

    def test_should_not_output_snp_with_non_canonical_character_in_alt(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATGGATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "..............K.........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)

    def test_should_not_output_snp_with_non_canonical_base_before_insertion(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATK*ATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "...............T........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)

    def test_should_not_output_snp_with_non_canonical_base_before_deletion(self):
        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            "ATCGATCGATCGATKTATCGATCGATCGATCGATCGATCG"
        ).with_read(
            "...............*........................", n_fwd=10, n_rev=10
        )

        expect = svc_driver.call()
        expect.with_output_vcf().record_count(0)
