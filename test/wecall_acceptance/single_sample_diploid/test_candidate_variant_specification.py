# All content Copyright (C) 2018 Genomics plc
from os import path
from unittest import expectedFailure

from wecall.bamutils.sample_bank import SampleBank
from wecall.bamutils.sequence_bank import AsciiVariantGenerator
from wecall.genomics.variant import Variant
from wecall.vcfutils.info_data import InfoData
from wecall.vcfutils.vcf_builder import VCFBuilder
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver
from wecall_test_drivers.variant_caller_builder import VariantCallerBuilderFromSampleBank
from wecall_test_drivers.variant_caller_wrapper import CANDIDATE_VARIANTS_FILE_KEY


class TestCandidateVariantSpecification(BaseTest):
    def calls_variants(self, ref, sequence_list, candidate_ascii_haplotypes, prior, expected_ascii_haplotypes):
        sample_bank = SampleBank(ref)
        sample_bank.add_sample_with_seqs_and_quals("TEST", sequence_list, 1, 0)

        variant_generator = AsciiVariantGenerator(sample_bank.reference)
        candidate_variants = variant_generator.get_variants(candidate_ascii_haplotypes)
        expected_variants = variant_generator.get_variants(expected_ascii_haplotypes)

        candidate_variant_list = VCFBuilder(path.join(self.work_dir, "candiate_variants.vcf"))
        candidate_variant_list.schema.set_info_data('AF', 'A', 'Float', 'Allele Frequency')
        for var in candidate_variants:
            candidate_variant_list.with_record_from_variant(
                var, info=InfoData(candidate_variant_list.schema, {"AF": prior})
            )
        candidate_variant_list.build().index()

        vc_wrapper_builder = VariantCallerBuilderFromSampleBank(sample_bank, self.work_dir)
        vc_wrapper_builder.configuration[CANDIDATE_VARIANTS_FILE_KEY] = candidate_variant_list.compressed_filename
        callset = vc_wrapper_builder.build().run().get_variant_callset(self)

        self.assertEqual(callset.get_variants(), set(expected_variants))

    @expectedFailure
    def test_incorrectly_aligned_reads_support_deletion(self):
        self.calls_variants(
            "TGTTATTAATCCCTTGTCAGATGTTATTAATCCCTTGTCAGT***CGCAAATATTTT",  # reference
            ["    ......................................GCAA           ",
             "    ......................................GCAA           ",
             "    ......................................GCAA           ",
             "    ......................................GCAA           ",
             "    ......................................GCAA           ",
             "    ......................................GCAA           ",
             ],
            # 01234567890123456789012345678901234567890123456789
            ["..........................................****..........."],  # candidates
            [0.72],
            ["..........................................****...........",  # expected output
             "..........................................****..........."]
        )

    def test_incorrectly_aligned_reads_supporting_snp(self):
        self.calls_variants(
            "TGTTATTAATCCCTTGCCAGAAACCATATCTTTTTTTTTTGCAAATATTTT",  # reference
            ["   TGTTATTAATCCCTTGCCATAAACCATATC                  ",  # input reads
             "   TAATCCCTTGCCATAAACCATATCTTTTTTTTTTG             "],

            # 01234567890123456789012345678901234567890123456789
            ["...................T..............................."],  # candidates
            [0.00000005],
            ["...................T...............................",  # expected output
             "...................T..............................."]
        )

    def test_incorrectly_aligned_reads_supporting_large_deletion_with_homopolymer_region(self):
        self.calls_variants(
            "TGTTATAAAAAAAAAAAATAAATAAATATATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCTCTAAAAAAAAT",
            ["...................................TCTCT........T  .T............T...T...T.T.T..................",
             # input reads that ALL support a long deletion
             "...................................T                                .T...T.T.T................. ",
             ],
            # 01234567890123456789012345678901234567890123456789
            ["...............................***********************************************.................."],
            # candidate
            [0.00064],
            ["...............................***********************************************..................",
             # expected output
             "...............................***********************************************.................."]
        )

    @expectedFailure  # "WIP Need to be able to modify the quality mapping"
    def test_should_not_call_snp_only_supported_by_low_quality_reads(self):
        self.calls_variants(
            "TGTTATTAATCCCTTGCCAGAAACCATATCTTTTTTTTTTGCAAATATTTT",  # reference
            ["........................T..........................",
             "                        1                          ",
             "........................T..........................",
             "                        1                          ",
             "........................T..........................",
             "                        1                          ",
             "........................T..........................",
             "                        1                          "],
            ["........................T.........................."],
            [1e-3],
            ["...................................................",
             "..................................................."]
        )

    def test_should_have_zero_bad_reads_for_candidate_variant_with_no_reads_covering_variant(self):
        chrom = "1"
        candidate_variant_list = VCFBuilder(path.join(self.work_dir, "candiate_variants.vcf"))
        candidate_variant_list.schema.set_info_data('AF', 'A', 'Float', 'Allele Frequency')
        variant_1 = Variant(chrom, 30, 'T', 'C')
        candidate_variant_list.with_record_from_variant(
            variant_1, info=InfoData(candidate_variant_list.schema, {"AF": [0.72]}))
        candidate_variant_list.build().index()

        svc_driver = SVCDriver(self)\
            .with_allow_MNP_calls(True)\
            .with_ref_sequence(
                "TGTTATTAATCCCTTGTCAGATGTTATTAATCCCTTGTCAGTCCCTTGTCAGT", chrom=chrom)\
            .with_read(
                "...........................C.. ......................", n_fwd=10, n_rev=10, sample_name='sample_1')\
            .with_read(
                "                                                     ", n_fwd=10, n_rev=10, sample_name='sample_2')\
            .with_candidate_variants_file(candidate_variant_list.compressed_filename)

        expect = svc_driver.call()

        vcf_expect = expect.with_output_vcf()
        vcf_expect.missing_record_for_variant(variant_1)
