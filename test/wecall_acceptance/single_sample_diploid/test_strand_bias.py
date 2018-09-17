# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestStrandBiasPValue(BaseTest):
    def test_should_get_unknown_value_if_all_reads_are_forward(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT", chrom="1"
        ).with_read(
            ".....................T....................", chrom="1", n_fwd=10, n_rev=0)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_record_for_variant(Variant("1", 21, "A", "T")) \
            .with_info().with_field("SBPV", [None])

    def test_should_get_unknown_value_if_all_reads_are_reverse(self):
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT", chrom="1"
        ).with_read(
            ".....................T....................", chrom="1", n_fwd=0, n_rev=10)

        vcf_expect = driver.call().with_output_vcf()
        vcf_expect.has_record_for_variant(Variant("1", 21, "A", "T")) \
            .with_info().with_field("SBPV", [None])


class TestCallingWithForwardAndReverseReads(AsciiWecallRunnerTest):

    def test_calls_snp_on_full_length_forward_reads(self):
        self.calls_variants(
            "ACGCCCCCTGCAAAAAAAAAA",
            ["...........C.........",
             ",,,,,,,,,,,c,,,,,,,,,",
             ",,,,,,,,,,,c,,,,,,,,,",
             "...........C........."],

            ["...........C.........",
             "...........C........."],  # Expected genotype
        )

    def test_calls_snp_on_long_forward_reads(self):
        self.calls_variants(
            "AAAAAAAAAAACGCCCCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ...................C.........           "]
        )

    def test_calls_snp_on_forward_and_reverse_reads(self):
        self.calls_variants(
            "AAAAAAAAAAACGCCCCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,           ",
             "    ,,,,,,,,,,t,,,,,,,,,,,,,,,,,,,,,,     ",
             "..............T.............              "]
        )

    def test_calls_del_and_snp_on_forward_and_reverse_reads(self):
        self.calls_variants(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,           ",
             "    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,     ",
             "..............*.............              "]
        )
