# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall.utils.interval import ChromInterval
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestRefCalling(BaseTest):
    def test_dont_call_reference_between_variant_and_insertion_due_to_vcf_rep_issues(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACG*CCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..........T................        ", chrom=chrom
        ).with_read(
            "  ...............T...........             ", chrom=chrom
        ).with_read(
            "    ...........T.*................        ", chrom=chrom
        ).with_read(
            "    ...........T.*.....................   ", chrom=chrom
        ).with_read(
            "...............T.*.........               ", chrom=chrom
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()

        # Has only 4 records which are:-
        vcf_expect.has_reference_calls_for_region(chrom, 0, 15)
        vcf_expect.has_record(chrom, 15, "C", "T")
        vcf_expect.has_record(chrom, 16, "G", "GT")
        vcf_expect.has_reference_calls_for_region(chrom, 17, 41)

    def test_calls_whole_region_as_reference_if_no_variants(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ...........................       ", chrom=chrom
        ).with_read(
            "  ...........................            ", chrom=chrom
        ).with_read(
            "    ...............................      ", chrom=chrom
        ).with_read(
            "    ...................................  ", chrom=chrom
        ).with_read(
            "...........................              ", chrom=chrom
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 41)

    def test_calls_correct_ref_calls_with_one_snp(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..........T................       ", chrom=chrom
        ).with_read(
            "  ...............T...........            ", chrom=chrom
        ).with_read(
            "    ...............................      ", chrom=chrom
        ).with_read(
            "    ...................................  ", chrom=chrom
        ).with_read(
            "...........................              ", chrom=chrom
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 17)
        vcf_expect.has_record(chrom, 17, "C", "T")
        vcf_expect.has_reference_calls_for_region(chrom, 18, 41)

    def test_calls_correct_ref_calls_with_one_mnp(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..........TTTT.............       ", chrom=chrom
        ).with_read(
            "  ...............TTTT........            ", chrom=chrom
        ).with_read(
            "    ...............................      ", chrom=chrom
        ).with_read(
            "    ...................................  ", chrom=chrom
        ).with_read(
            "...........................              ", chrom=chrom
        ).with_output_ref_calls(True).with_allow_MNP_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 17)
        vcf_expect.has_record(chrom, 17, "CCCC", "TTTT")
        vcf_expect.has_reference_calls_for_region(chrom, 21, 41)

    def test_calls_ref_calls_correctly_with_mnp_that_contains_snp(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..........TTTT.............       ", chrom=chrom
        ).with_read(
            "  ...............TTTT........            ", chrom=chrom
        ).with_read(
            "    ..............T................      ", chrom=chrom
        ).with_read(
            "    ..............T....................  ", chrom=chrom
        ).with_read(
            "..................T........              ", chrom=chrom
        ).with_output_ref_calls(True).with_allow_MNP_calls(True)

        vcf_expect = driver.call().with_output_vcf()

        vcf_expect.has_reference_calls_for_region(chrom, 0, 17)
        vcf_expect.has_record(chrom, 17, "CCCC", "TTTT")
        vcf_expect.has_record(chrom, 18, "C", "T")
        vcf_expect.has_reference_calls_for_region(chrom, 21, 41)

    def test_calls_correct_ref_calls_with_one_del(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACTCCCCATAAAAAAATTTTTTTTTTT",
        ).with_read(
            "       ..........*................       ", chrom=chrom
        ).with_read(
            "  ...............*...........            ", chrom=chrom
        ).with_read(
            "    ...............................      ", chrom=chrom
        ).with_read(
            "    ...................................  ", chrom=chrom
        ).with_read(
            "...........................              ", chrom=chrom
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 16)
        vcf_expect.has_record(chrom, 16, "TC", "T")
        vcf_expect.has_reference_calls_for_region(chrom, 18, 41)

    def test_calls_correct_ref_calls_with_one_ins(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACG*CCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..........*................        ", chrom=chrom
        ).with_read(
            "  ...............*...........             ", chrom=chrom
        ).with_read(
            "    .............T.................       ", chrom=chrom
        ).with_read(
            "    .............T.....................   ", chrom=chrom
        ).with_read(
            ".................T.........               ", chrom=chrom
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 16)
        vcf_expect.has_record(chrom, 16, "G", "GT")
        vcf_expect.has_reference_calls_for_region(chrom, 17, 41)

    def test_calls_correct_ref_calls_with_cluster_of_variants(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAATAACGCACG*CCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..*.......*................        ", chrom=chrom
        ).with_read(
            "  .......*.......*...........             ", chrom=chrom
        ).with_read(
            "    .............T......T..........       ", chrom=chrom
        ).with_read(
            "    .............T......T..............   ", chrom=chrom
        ).with_read(
            ".................T......T..               ", chrom=chrom
        ).with_output_ref_calls(True)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_reference_calls_for_region(chrom, 0, 8)
        vcf_expect.has_record(chrom, 8, "TA", "T")
        vcf_expect.has_reference_calls_for_region(chrom, 10, 16)
        vcf_expect.has_record(chrom, 16, "G", "GT")
        vcf_expect.has_reference_calls_for_region(chrom, 17, 23)
        vcf_expect.has_record(chrom, 23, "A", "T")
        vcf_expect.has_reference_calls_for_region(chrom, 24, 41)

    def test_calls_reference_on_location_with_low_quality_variant_support(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAATAACGCACGCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..*.......................        ", n_fwd=2, n_rev=1
        ).with_read(
            ".................T.....T.........        ",
            "                 1                       ", n_fwd=1, n_rev=1
        ).with_read(
            ".......................T..               ", n_fwd=1, n_rev=0
        ).with_output_ref_calls(True)

        expect = driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.has_reference_calls(ChromInterval(chrom, 0, 8))

        vcf_expect.has_record_for_variant(Variant(chrom, 8, "TA", "T"))

        vcf_expect.has_reference_calls(ChromInterval(chrom, 10, 23))

        vcf_expect.has_record_for_variant(Variant(chrom, 23, "A", "T"))

        vcf_expect.has_reference_calls(ChromInterval(chrom, 24, 41))

    def test_calls_correct_reference_between_clusters_with_uncalled_indel_between(self):
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAATAACGCACGCCCCATAAAAAAATTTTTTTTTTT", chrom=chrom
        ).with_read(
            "       ..*.......................        ", n_fwd=10, n_rev=10
        ).with_read(
            ".......................*............     ", n_fwd=1, n_rev=1
        ).with_read(
            ".......................T..               ", n_fwd=10, n_rev=10
        ).with_output_ref_calls(True).with_max_cluster_distance(5)

        expect = driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect.has_reference_calls(ChromInterval(chrom, 0, 8))

        vcf_expect.has_record_for_variant(Variant(chrom, 8, "TA", "T"))

        vcf_expect.has_reference_calls(ChromInterval(chrom, 10, 23))

        vcf_expect.has_record_for_variant(Variant(chrom, 23, "A", "T"))

        vcf_expect.has_reference_calls(ChromInterval(chrom, 24, 41))
