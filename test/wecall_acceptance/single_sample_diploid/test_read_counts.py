# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall.genomics.variant import Variant
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestVariantReadCounts(AsciiWecallRunnerTest):
    def test_should_count_reads_that_do_not_overlap_the_calling_region(self):
        sample = "seed"
        chrom = "1"

        driver = SVCDriver(self).with_allow_MNP_calls(True)
        driver.with_ref_sequence(
            "GAAAAAAAAAAACGCACCCCCAAATTTTTTTTAA***********AAAATAAAAAACGCACCCCCAAATTTTTTTTAA", chrom=chrom
        ).with_read(
            "                                             ..........G......................",
            n_fwd=10, n_rev=10, sample_name=sample,
        ).with_read(
            "..................................AAAATAAAAAG                                 ",
            n_fwd=10, n_rev=10, sample_name=sample,
        ).with_region_string("{}:{}-{}".format(chrom, 34, 55))

        vcf = driver.call().with_output_vcf()
        vcf \
            .has_record_for_variant(Variant(chrom, 44, 'A', 'G')) \
            .with_info() \
            .with_field("DP", [40]) \
            .with_field("VC", [40])

    def test_should_record_the_read_support_mnp_and_snp(self):
        sample = "seed"
        chrom = "1"

        driver = SVCDriver(self).with_allow_MNP_calls(True)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCAAATTTTTTTTAA", chrom=chrom
        ).with_read(
            "           ......................", n_fwd=0, n_rev=1, sample_name=sample,
        ).with_read(
            "..........G......................", n_fwd=10, n_rev=10, sample_name=sample,
        ).with_read(
            "..........CT.....................", n_fwd=10, n_rev=9, sample_name=sample)

        vcf = driver.call().with_output_vcf()

        vcf \
            .with_samples([sample]) \
            .record_count(2)

        vcf \
            .has_record_for_variant(Variant(chrom, 10, 'AC', 'CT')) \
            .with_sample(sample) \
            .has_read_depth(40) \
            .has_allelic_read_support(1, 19) \
            .has_variant_allelic_frequency(19.0 / 40.0)

        vcf \
            .has_record_for_variant(Variant(chrom, 10, 'A', 'G')) \
            .with_sample(sample) \
            .has_read_depth(39) \
            .has_allelic_read_support(0, 20) \
            .has_variant_allelic_frequency(20.0 / 39.0)

    def test_should_record_the_read_support_insertion_and_snp_on_same_strand(self):
        sample = "sample"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAACGTAGCTG*GCACCCCCAAA", chrom=chrom
        ).with_read(
            "...........*...........", n_fwd=10, n_rev=10, sample_name=sample,
        ).with_read(
            "..........CT...........", n_fwd=10, n_rev=10, sample_name=sample,
        )

        vcf = driver.call().with_output_vcf()

        vcf \
            .with_samples([sample]) \
            .record_count(2)

        vcf \
            .has_record_for_variant(Variant(chrom, 9, 'T', 'TC')) \
            .with_sample(sample) \
            .has_read_depth(40) \
            .has_allelic_read_support(20, 20) \
            .has_genotype("0/1")

    def test_should_record_the_read_support_deletion_and_snp_on_same_strand(self):
        sample = "sample"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAACGTAGCTGTGCACCCCCAAA", chrom=chrom
        ).with_read(
            ".......................", n_fwd=10, n_rev=10, sample_name=sample,
        ).with_read(
            "..........C*...........", n_fwd=10, n_rev=10, sample_name=sample,
        )

        vcf = driver.call().with_output_vcf()

        vcf \
            .with_samples([sample]) \
            .record_count(2)

        vcf \
            .has_record(chrom, 9, 'TG', 'T') \
            .with_sample(sample) \
            .has_read_depth(40) \
            .has_allelic_read_support(20, 20) \
            .has_genotype("0/1")

    def test_should_report_unkown_value_for_allele_frequence_when_depth_is_zero(self):
        # Only way to output depth zero for sample and variant is to have
        # another good sample
        good_sample = "good_sample"
        empty_sample = "empty_sample"
        chrom = "1"

        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAACGTAGCTGTGCACCCCCAAA", chrom=chrom
        ).with_read(
            "..........T............", n_fwd=10, n_rev=10, sample_name=good_sample,
        ).with_read(
            ".......................", n_fwd=0, n_rev=0, sample_name=empty_sample,
        )
        vcf = driver.call().with_output_vcf()
        vcf \
            .with_samples([good_sample, empty_sample]) \
            .record_count(1)

        vcf \
            .has_record_for_variant(Variant(chrom, 10, "G", "T")) \
            .with_sample(empty_sample) \
            .has_read_depth(0) \
            .has_variant_allelic_frequency(None)

    def test_single_snp_in_clean_data(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ...................C........            ",
             "    .................C..............      ",
             "    .................C..................  ",
             ".....................C......              "],
            [(21, "A", "C", {"DP": [5], "AD": [0, 5]})]
        )

    def test_single_snp_in_dirty_data(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............T.............       ",
             "                     2                    ",
             "  ...................T........            ",
             "                     2                    ",
             "    .................C..............      ",
             "    .................C..................  ",
             ".....................C......              "],
            [(21, "A", "C", {"DP": [5], "AD": [0, 3]})], config_dict={"minBaseQual": 20}
        )

    def test_snps_at_same_location_should_have_expected_read_counts(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ...................C........            ",
             "    .................C..............      ",
             "    .................C..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "C", {"DP": [8], "AD": [0, 4]}),
             (21, "A", "T", {"DP": [8], "AD": [0, 4]})]
        )

    def test_snp_and_mnp_at_same_location_with_mnp_containing_snp_should_have_expected_read_counts(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............TC............       ",
             "  ...................TC.......            ",
             "    .................TC.............      ",
             "    .................TC.................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "AT", "TC", {"DP": [8], "AD": [0, 4]}),
             (21, "A", "T", {"DP": [8], "AD": [0, 4]})], config_dict={"allowMNPCalls": "True"}
        )

    def test_should_only_count_reads_that_overlap_the_snp(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTTAT",
            [" .....................                      ",
             "..................... ..................... ",
             "       ..............T.............         ",
             "  ...................T........              ",
             "    .................T..............        ",
             "    .................T..................    ",
             "    .................T..................    ",
             "    .................T..................    ",
             "    .................T..................    ",
             ".....................T......                "],
            [(21, "A", "T", {"DP": [9], "AD": [1, 8]})]
        )

    def test_should_only_count_reads_that_overlap_the_ins(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCC*CAAAAAAAATTTTTTTTTTTAT",
            [" ....................*.                     ",
             "..................... ..................... ",
             "       ..............T.............         ",
             "  ...................T........              ",
             "    .................T..............        ",
             "    .................T..................    ",
             "    .................T..................    ",
             "    .................T..................    ",
             "    .................T..................    ",
             ".....................T......                "],
            [(20, "C", "CT", {"DP": [9], "AD": [1, 8]})]
        )

    def test_should_only_count_reads_that_overlap_the_del(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTTAT",
            [" .....................                      ",
             "..................... ..................... ",
             "       ..............*.............         ",
             "  ...................*........              ",
             "    .................*..............        ",
             "    .................*..................    ",
             "    .................*..................    ",
             "    .................*..................    ",
             "    .................*..................    ",
             ".....................*......                "],
            [(20, "CA", "C", {"DP": [9], "AD": [1, 8]})]
        )

    def test_snps_at_same_location_should_have_expected_read_counts_with_anonymous_reads(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ...................C........            ",
             "  ............................            ",
             "    .................C..............      ",
             "    .................C..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             "    .................T..................  ",
             ".....................T......              "],
            [(21, "A", "C", {"DP": [9], "AD": [1, 4]}),
             (21, "A", "T", {"DP": [9], "AD": [1, 4]})]
        )

    def test_should_have_expected_sample_info_data_when_half_reads_support_snp(self):
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ............................            ",
             "    .................C..............      ",
             "    ....................................  ",
             "    ....................................  ",
             "    .................C..................  ",
             "    ....................................  ",
             "    .................C..................  ",
             "    ....................................  ",
             ".....................C......              "],
            [(21, "A", "C", {"DP": [10], "AD": [5, 5]})]
        )

    def test_all_reads_support_deletion(self):
        self.calls_variants_with_sample_data(
            "TGTTATTAATCCCTTGTCAGATTTTTTTTTTGCAAATATTTT",
            ["       ..............**********........   ",
             "  ...................**********........   ",
             "    .................**********.....      ",
             "    .................**********.........  ",
             "    .................**********.........  ",
             "    .................**********.........  ",
             "    .................**********.........  ",
             "    .................**********.........  ",
             "    .................**********.........  ",
             ".....................**********.......... "],
            [(20, "ATTTTTTTTTT", "A", {"DP": [10], "AD": [0, 10]})]
        )

    def test_incorrectly_aligned_reads_support_deletion(self):
        self.calls_variants_with_sample_data(
            "TGTTATTAATCCCTTGTCAGA*********TTTTTTTTTTGCAAATATTTTCTGATGAGTACGG",
            ["       ..............GCAAATATT                                  ",
             "  ...................GCAAATATT                                  ",
             "    .................GCAAATATT                                  ",
             "    .................*******************...................     ",
             "    .................*******************...................     ",
             "    .................*******************...................     ",
             "    .................*******************...................     ",
             "    .................*******************...................     ",
             "    .................*******************...................     ",
             ".....................*******************..............          "],
            [(20, "ATTTTTTTTTT", "A", {"DP": [10], "AD": [0, 10]})]
        )

    def test_should_only_count_reads_that_properly_overlap_a_deletion_as_support(self):
        self.calls_variants_with_sample_data(
            "TTGTATTTCCTGTTATTAATCCCTTGTCAGATTTTTTTTTTGCAAATATTTT",
            ["...............................**********....       ",
             "...............................**********....       ",
             "...............................**********...        ",
             "...............................**********..         ",
             "...............................**********.          ",
             "...............................**********           "],
            [(30, "ATTTTTTTTTT", "A", {"DP": [6], "AD": [1, 5]})]
        )

    def test_should_count_reads_that_partially_overlap_an_insertion_if_inserted_sequence_is_distinct_from_subsequent_reference_sequence(self):  # noqa
        self.calls_variants_with_sample_data(
            "AAAAAAAAAAACGCACCCCCC****ATAAAAAAAATTTTTTTTTTT",
            ["       ..............CTCT.............        ",
             ".....................CTC                      ",
             "    .................CTCT..............       ",
             "    .................CTCT...............      ",
             ".....................CTC                      ",
             "    .................CTCT..................   ",
             ".....................CTC                      ",
             "    .................CTCT..................   ",
             ".....................CTC                      ",
             ".....................CTCT......               "],
            [(20, "C", "CCTCT", {"DP": [10], "AD": [0, 10]})]
        )

    def test_should_count_reads_support_consistent_with_left_alignment(self):
        self.calls_variants_with_sample_data(
            "TTAGATTTTTTTCCTATAGAATTG**TTTTTTTTTTTT**ATTTCCTGTTATTAATCCCTTGTCAGATTTTTT",
            ["........................TT............**.................................",
             "........................**............TT.................................",
             "........................**............TT.................................",
             ],
            [(23, "G", "GTT", {"DP": [3], "AD": [0, 3]})]
        )

    @expectedFailure
    def test_should_count_read_support_consistent_with_large_repeat_motif_across_different_reads(self):
        self.calls_variants_with_sample_data(
            "GAATTGTTTTTTTTTTTTATTTCCTGTTTTTTTTTTTTATTTCCTGTT*********A*********TTCTGTC*********CATTCTGTCCGGATTCAGATACTTTGCCCATTTTTAAGTTGGATCATTAGATTTTTTTCCTAT",  # noqa
            ["................................................ATTCTGTCC.*********.......*********...                                                            ",  # noqa
             "................................................ATTCTGTCC.*********.......*********...                                                            ",  # noqa
             "    ............................................*********.TTCTGTCCA.......*********.......                                                        ",  # noqa
             "    ............................................*********.TTCTGTCCA.......*********.......                                                        ",  # noqa
             "                                                         .*********.......CATTCTGTC...............................................................",  # noqa
             "                                                         .*********.......CATTCTGTC..............................................................."  # noqa
             ],
            [(47, "T", "TATTCTGTCC", {"DP": [6], "AD": [0, 6]})]
        )

    def test_should_count_partial_support_for_non_completely_overlapping_reads(self):
        self.calls_variants_with_sample_data(
            "GAATTGTTTTTTTTTTTTATTTCCTGTTTTTTTTTTTTATTTTCCTGTTATTCTGTCCATTCTGTCCGGATTCAGATACTTTGCCCATTT",
            [".......................................A..                                                ",
             "                        ...............A..G..........................                     ",
             "                        ...............A..G..........................                     ",
             "                                        ..G.............................................. "],
            [(39, "TTTT", "ATTG", {"DP": [4], "AD": [0, 3]})], config_dict={"allowMNPCalls": "True"}
        )

    def test_should_count_reads_that_partially_overlap_insertions_as_supporting_reads(self):
        self.calls_variants_with_sample_data(
            "GAATTGTTTTTTTTTTTTATTTCCTGTTTTTTTTTTTTAT**********ATTCCTGTTATTCTGTCCATTCTGTCCGGATTCAGATACTTTGCCCAT",
            ["........................................TATTCTGTCC..                                              ",
             "........................................TATTCTGTCC..                                              ",
             "........................................TATTCTGT                                                  ",
             "........................................TATTCTGT                                                  ",
             "........................................TATTCTG                                                   ",
             "........................................TATTCTG                                                   ",
             "........................................TATTCT                                                    ",
             "........................................TATTCT                                                    ",
             "........................................TATTC                                                     ",
             "........................................TATTC                                                     ",
             "........................................TATT                                                      ",
             "........................................TATT                                                      ",
             "........................................TAT                                                       ",
             "........................................TAT                                                       ",
             "........................................TA                                                        ",
             "........................................TA                                                        ",
             "........................................T                                                         ",
             "........................................T                                                         "],
            [(39, "T", "TTATTCTGTCC", {"DP": [18], "AD": [0, 18]})]
        )

    def test_should_count_reads_that_support_both_alternate_and_reference_allele_as_half(self):
        self.calls_variants_with_sample_data(
            "GAATTGTTTTTTTTTTTTATTTCCTGTTTTTTTTTTTTAT**********TTTCCTGTTATTCTGTCCATTCTGTCCGGATTCAGATACTTTGCCCAT",
            ["........................................TATTCTGTCC..                                              ",
             "........................................TATTCTGTCC..                                              ",
             "........................................T                                                         ",
             "........................................T                                                         "],
            [(39, "T", "TTATTCTGTCC", {"DP": [4], "AD": [1, 3]})]
        )


class TestReadCountsInVariantCluster(AsciiWecallRunnerTest):

    def test_coverage_and_supporting_counts_for_one_snp(self):
        self.calls_variants_with_sample_data(
            "TGCGAATACATCGCACCCCCCATACAACAAATTTGTCTATTG",
            ["       ..............C.............       ",
             "  ............................            ",
             "    .................C..............      ",
             "    ....................................  ",
             "    ....................................  ",
             "    .................C..................  ",
             "    ....................................  ",
             "    .................C..................  ",
             "    ....................................  ",
             ".....................C......              "],
            [(21, "A", "C", {"DP": [10], "AD": [5, 5]})]
        )

    def test_coverage_and_supporting_counts_for_two_snps(self):
        self.calls_variants_with_sample_data(
            "TGCGAATACATCGCACCCCCCATACAACAAATTTGTCTATTG",
            ["       ..............C.............       ",
             "  ..............A.............            ",
             "    .................C..............      ",
             "    ............A.......................  ",
             "    ............A.......................  ",
             "    .................C..................  ",
             "    ............A.......................  ",
             "    .................C..................  ",
             "    ............A.......................  ",
             ".....................C......              "],
            [(16, "C", "A", {"DP": [10], "AD": [5, 5]}),
             (21, "A", "C", {"DP": [10], "AD": [5, 5]})]
        )

    def test_coverage_and_supporting_counts_for_two_snps_in_one_cluster_when_some_reads_only_cover_one_snp_position(self):  # noqa
        self.calls_variants_with_sample_data(
            "ACCTGACCTGCGAATACATCGCACCCCCCATCTGCTACAACAAATTTGTCTATTGCGTAATGCC",
            [".............................C..     A..........................",
             "                                     A..........................",
             "                                     A..........................",
             "                                     A..........................",
             "                                     A..........................",
             ".............................C..                                ",
             ".............................C..                                ",
             ".............................C..                                ",
             ".............................C..                                "],
            [(29, "A", "C", {"DP": [5], "AD": [0, 5]}),
             (37, "C", "A", {"DP": [5], "AD": [0, 5]})], config_dict={"allowMNPCalls": False}
        )

    def test_coverage_and_supporting_counts_for_two_snps_in_two_clusters_when_some_reads_only_cover_one_snp_position(self):  # noqa
        self.calls_variants_with_sample_data(
            "ACCTGACCTGCGAATACATCGCACCCCCCATCCTGCCTGCTCTCCAACAAATTTGTCTATTGCGTAATGCC",
            ["                                            A..........................",
             "                                            A..........................",
             "                                            A..........................",
             "                                            A..........................",
             "                                            A..........................",
             ".............................C..                                       ",
             ".............................C..                                       ",
             ".............................C..                                       ",
             ".............................C..                                       ",
             ".............................C..                                       "],
            [(29, "A", "C", {"DP": [5], "AD": [0, 5]}),
             (44, "C", "A", {"DP": [5], "AD": [0, 5]})], config_dict={"allowMNPCalls": False}
        )
