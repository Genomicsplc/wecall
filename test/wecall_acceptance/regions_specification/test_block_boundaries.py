# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest


class TestBlockBoundariesCluster(AsciiWecallRunnerTest):
    def calls_variants_with_defined_block_size(self, ref, sequence_list, expected_variants=None):
        self.calls_variants(
            ref,
            sequence_list,
            config_dict={
                "maxBlockSize": self.max_block_size,
                "maxClusterDist": 3,
                "minClusterDist": 3},
            expected_ascii_haplotypes=None,
            expected_variant_stubs=expected_variants,
            n_fwd=20,
            n_rev=20)

    def test_should_call_del_and_snp_either_side_of_border(self):
        self.max_block_size = 20
        self.calls_variants_with_defined_block_size(
            # 0123456789012345678901234567890123456789"
            "ACGCTCACTGACGCTCACTGATACTGACTGATCGCTGGTT",
            ["T................*G*T*C................."],  # input

            [
                (0, "A", "T"), (16, "ACT", "A"), (19, "GA", "G"), (22, "A", "C"),
            ]
        )

    def test_should_call_snp_at_the_end_of_the_reference(self):
        self.max_block_size = 10
        self.calls_variants_with_defined_block_size(
            # 012345678901234567890123456789"
            "ACGCTCACTGATACTGACTGATCGCTGGTT",
            [".............................C"],  # input

            [(29, "T", "C"), ]
        )

    def test_should_call_snp_at_the_begin_of_the_reference(self):
        self.max_block_size = 10
        self.calls_variants_with_defined_block_size(
            # 012345678901234567890123456789"
            "ACGCTCACTGATACTGACTGATCGCTGGTT",
            ["T............................."],  # input

            [(0, "A", "T"), ]
        )

    def test_should_left_align_isolated_del_across_boundary_when_there_are_nearby_variants_on_left(self):
        self.max_block_size = 20
        self.calls_variants_with_defined_block_size(
            # 0123456789012345678901234567890123456789"
            "ACGCTCACTTACGCTCACTTTTTCTGACTGATCGCTGGTT",
            ["T...............T.....*................."],  # input

            [(0, "A", "T"), (16, "A", "T"), (17, "CT", "C")]
        )
