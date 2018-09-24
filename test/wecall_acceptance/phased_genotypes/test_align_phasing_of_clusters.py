# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver

MAX_PHRED = 3000


class TestAlignPhasingOfClusters(BaseTest):
    def test_phase_alignment_for_two_snps_in_different_clusters_on_different_strands(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                ".......A.................................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".............................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_max_cluster_distance(10) \

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect.has_record_for_variants(Variant(chrom, 7, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("7")

        vcf_expect.has_record_for_variants(Variant(chrom, 29, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("7")

    def test_phase_alignment_for_two_snps_in_different_clusters_on_different_strands_for_reference_calling(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                ".......A.................................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".............................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_max_cluster_distance(20) \
            .with_output_ref_calls(1) \
            .with_verbosity(6)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(7)

        vcf_expect.has_record_for_variants(Variant(chrom, 7, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("7")

        vcf_expect.has_record_for_variants(Variant(chrom, 29, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("7")

    def test_phase_alignment_for_two_snps_in_different_clusters_on_same_strands(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                ".......A.....................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".........................................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(20)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect.has_record_for_variants(Variant(chrom, 7, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("7")

        vcf_expect.has_record_for_variants(Variant(chrom, 29, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("7")

    def test_phase_not_aligns_for_hom_snp_in_first_cluster(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                ".......A.....................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".......A.................................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(10)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect.has_record_for_variants(Variant(chrom, 7, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("7")

        vcf_expect.has_record_for_variants(Variant(chrom, 29, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("28")

    def test_phase_aligns_for_hom_snp_in_second_cluster(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                ".......A.....................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".............................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(10)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect.has_record_for_variants(Variant(chrom, 7, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("7")

        vcf_expect.has_record_for_variants(Variant(chrom, 29, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("7")

    def test_phase_aligns_for_hom_snp_in_both_cluster(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                ".......A.....................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                ".......A.....................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(10)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect.has_record_for_variants(Variant(chrom, 7, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("7")

        vcf_expect.has_record_for_variants(Variant(chrom, 29, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("7")

    def test_phase_alignment_for_hom_and_het_variants_on_different_clusters(self):
        sample_name = "a_sample"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                ".......A...A.................A...........", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "...........A...............A.............", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(10) \
            .with_min_cluster_distance(5)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(4)

        # phase set id is start of cluster
        vcf_expect.has_record_for_variants(Variant(chrom, 7, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("5")

        vcf_expect.has_record_for_variants(Variant(chrom, 11, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("5")

        vcf_expect.has_record_for_variants(Variant(chrom, 27, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("5")

        vcf_expect.has_record_for_variants(Variant(chrom, 29, "G", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("5")

    def test_phase_alignment_for_het_variants_for_two_samples(self):
        sample_name_1 = "sample1"
        sample_name_2 = "sample2"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                "...........A.............................", n_fwd=10, n_rev=10, sample_name=sample_name_1) \
            .with_read(
                "...........................A.............", n_fwd=10, n_rev=10, sample_name=sample_name_1) \
            .with_read(
                "...........A.............................", n_fwd=10, n_rev=10, sample_name=sample_name_2) \
            .with_read(
                "...........................A.............", n_fwd=10, n_rev=10, sample_name=sample_name_2) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(10) \
            .with_min_cluster_distance(5)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        # phase set id is start of cluster
        vcf_expect.has_record_for_variants(Variant(chrom, 11, "G", "A"))\
            .with_sample(sample_name_1)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("6")

        vcf_expect.has_record_for_variants(Variant(chrom, 27, "G", "A"))\
            .with_sample(sample_name_1)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("6")

        vcf_expect.has_record_for_variants(Variant(chrom, 11, "G", "A"))\
            .with_sample(sample_name_2)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("6")

        vcf_expect.has_record_for_variants(Variant(chrom, 27, "G", "A"))\
            .with_sample(sample_name_2)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("6")

    def test_phase_alignment_for_het_variants_for_three_clusters(self):
        sample_name = "sample1"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                "....A.............C......................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "..................................T......", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(5) \
            .with_min_cluster_distance(5)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(3)

        vcf_expect.has_record_for_variants(Variant(chrom, 4, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 18, "T", "C"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 34, "A", "T"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("4")

    def test_phase_alignment_for_het_variants_for_three_clusters_when_middle_cluster_is_homozygous(self):
        sample_name = "sample1"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                "....A.............C......................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "..................C...............T......", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(5) \
            .with_min_cluster_distance(5)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(3)

        vcf_expect.has_record_for_variants(Variant(chrom, 4, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 18, "T", "C"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 34, "A", "T"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("29")

    def test_phase_alignment_for_het_variants_for_three_clusters_when_first_cluster_is_homozygous(self):
        sample_name = "sample1"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                "....A....................................", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_read(
                "....A.............C...............T......", n_fwd=10, n_rev=10, sample_name=sample_name) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(5) \
            .with_min_cluster_distance(5)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(3)

        vcf_expect.has_record_for_variants(Variant(chrom, 4, "C", "A"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 18, "T", "C"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("13")

        vcf_expect.has_record_for_variants(Variant(chrom, 34, "A", "T"))\
            .with_sample(sample_name)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("13")

    def test_phase_alignment_for_het_variants_for_two_samples_with_missing_read_support(self):
        sample_name_1 = "sample1"
        sample_name_2 = "sample2"
        chrom = "1"

        svc_driver = SVCDriver(self)\
            .with_ref_sequence(
                "ACGCCCCCTGGGGGGGGGTGGGGGGGGGGGCAAAAAAAAAA", chrom=chrom) \
            .with_read(
                "...........A.............................", n_fwd=10, n_rev=10, sample_name=sample_name_1) \
            .with_read(
                "...........................A.............", n_fwd=10, n_rev=10, sample_name=sample_name_1) \
            .with_read(
                "...........A.............                ", n_fwd=10, n_rev=10, sample_name=sample_name_2) \
            .with_read(
                ".........................                ", n_fwd=10, n_rev=10, sample_name=sample_name_2) \
            .with_output_phased_genotypes(True) \
            .with_allow_MNP_calls(False) \
            .with_max_cluster_distance(5) \
            .with_min_cluster_distance(5)

        vcf_expect = svc_driver.call()\
            .with_output_vcf()\
            .record_count(2)

        vcf_expect.has_record_for_variants(Variant(chrom, 11, "G", "A"))\
            .with_sample(sample_name_1)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("6")

        vcf_expect.has_record_for_variants(Variant(chrom, 27, "G", "A"))\
            .with_sample(sample_name_1)\
            .has_exact_phased_genotypes("0|1")\
            .has_phase_set_id("6")

        vcf_expect.has_record_for_variants(Variant(chrom, 11, "G", "A"))\
            .with_sample(sample_name_2)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("6")

        vcf_expect.has_record_for_variants(Variant(chrom, 27, "G", "A"))\
            .with_sample(sample_name_2)\
            .has_exact_phased_genotypes("0|0")\
            .has_phase_set_id("6")

    def test_phase_alignment_for_overlapping_variants_for_two_samples(self):
        chrom = "1"
        sample1 = "sample1"
        sample2 = "sample2"

        svc_driver = SVCDriver(self)
        svc_driver\
            .with_ref_sequence(
                'CGAGCGATACAGATAAAGACATCGAGTGA', chrom=chrom) \
            .with_read(
                '....T................A.......', n_rev=5, n_fwd=0, chrom=chrom, sample_name=sample1) \
            .with_read(
                '.....................A.......', n_rev=5, n_fwd=0, chrom=chrom, sample_name=sample1) \
            .with_read(
                '....T................A.......', n_rev=5, n_fwd=0, chrom=chrom, sample_name=sample2) \
            .with_read(
                '.....................C.......', n_rev=5, n_fwd=0, chrom=chrom, sample_name=sample2)\
            .with_max_cluster_distance(8) \
            .with_min_cluster_distance(8) \
            .with_verbosity(6)

        vcf_expect = svc_driver.call() \
            .with_output_vcf() \
            .record_count(3)

        vcf_expect.has_record_for_variants(Variant(chrom, 4, "C", "T"))\
            .with_sample(sample1)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 21, "T", "C"))\
            .with_sample(sample1)\
            .has_exact_phased_genotypes(".|.")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 21, "T", "A"))\
            .with_sample(sample1)\
            .has_exact_phased_genotypes("1|1")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 4, "C", "T"))\
            .with_sample(sample2)\
            .has_exact_phased_genotypes("1|0")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 21, "T", "C"))\
            .with_sample(sample2)\
            .has_exact_phased_genotypes(".|1")\
            .has_phase_set_id("4")

        vcf_expect.has_record_for_variants(Variant(chrom, 21, "T", "A"))\
            .with_sample(sample2)\
            .has_exact_phased_genotypes("1|.")\
            .has_phase_set_id("4")
