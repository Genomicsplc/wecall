# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall.bedutils.bedwriter import BEDWriterContextManager, BEDIndexer
import os
from wecall.genomics.reference_chromosome import DEFAULT_CHROM
from wecall.utils.interval import Interval, ChromInterval
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest


class TestCallingInsideBedRegion(AsciiWecallRunnerTest):
    def calls_variants_in_bed_region(
            self,
            ref,
            sequence_list,
            expected_ascii_haplotypes=None,
            expected_variant_stubs=None,
            bed_regions=[],
            compress=False
    ):

        bed_filename = os.path.join(self.work_dir, "test.bed")
        with BEDWriterContextManager(bed_filename) as bed_writer:
            # TODO - link chromosomes
            for region in bed_regions:
                bed_writer.write_chrom_interval(ChromInterval(DEFAULT_CHROM, region.start, region.end))

        if compress:
            indexer = BEDIndexer(bed_filename)
            indexer.index()
            bed_filename = indexer.compressed_filename

        self.calls_variants(
            ref,
            sequence_list,
            expected_ascii_haplotypes,
            expected_variant_stubs,
            n_fwd=10,
            n_rev=10,
            config_dict={
                "regions": bed_filename,
                "allowMNPCalls": "True"})

    def test_should_call_variants_when_whole_read_within_region(self):
        self.calls_variants_in_bed_region(
            "ACGCCCCCTGCAAAAAAAAAA",
            [".T...................",
             "...........C........."],
            [".T...................",
             "...........C........."],

            bed_regions=[Interval(0, 12)]
        )

    def test_should_call_variants_when_whole_read_within_region_with_compressed_bed_file(self):
        self.calls_variants_in_bed_region(
            "ACGCCCCCTGCAAAAAAAAAA",
            [".T...................",
             "...........C........."],
            [".T...................",
             "...........C........."],

            bed_regions=[Interval(0, 12)], compress=True
        )

    @expectedFailure
    def test_should_call_only_variant_that_is_within_region(self):
        self.calls_variants_in_bed_region(
            "AAAAAAAAACGCCCCCTGCAAAAAAAAAA",
            [".........T...................",
             "..................T.........."],

            [".........T...................",
             ".............................", ],

            bed_regions=[Interval(8, 15)],
        )

    def test_should_call_variants_with_touching_and_unsorted_regions(self):
        self.calls_variants_in_bed_region(
            "ACGCCCCCTGCAAAAAAAAAA",
            ["......T...T.........."],
            expected_variant_stubs=[(6, "CCTGC", "TCTGT")],
            bed_regions=[Interval(7, 12), Interval(2, 7), Interval(4, 5)]
        )
