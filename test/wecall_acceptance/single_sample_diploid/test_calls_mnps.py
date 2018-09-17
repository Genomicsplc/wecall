# All content Copyright (C) 2018 Genomics plc
from unittest import expectedFailure

from wecall.genomics.variant import Variant
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestMNPCalling(AsciiWecallRunnerTest):
    @expectedFailure
    def test_should_call_mnp_at_the_end_of_read(self):
        expected_variant_stubs = {
            (3, "CTT", "C"),
            (18, "GT", "G"),
            (25, "C", "T"),
            # this MNP is not called and neither are the composing SNPs
            (37, "GCCG", "ACCT"),
            (45, "CTT", "C"),
            (53, "T", "TCTG")
        }

        self.calls_variants(
            "AACCTTGGACGTTATTCTGTCAATGCATCCCATTGCCGCCGCAACCTTGGACGT***TATTCTGTC",
            [" ...**.............*.....T...........A..T. ...**......CTG.........", ],
            n_fwd=10, n_rev=10,
            expected_variant_stubs=expected_variant_stubs
        )

    def test_calls_mnp_formed_by_overlapping_reads(self):
        sn = "a_sample"

        svc_driver = SVCDriver(self).with_allow_MNP_calls(True)
        svc_driver.with_ref_sequence(
            "AACCTTGGACGTTATTCTGTCAATGCATCCCATTGCCGCCGCAACCTTGGACGTTATTCTGTC", chrom="1"
        ).with_read(
            "..................T.. ....C.......C............................", sample_name=sn, n_fwd=3, n_rev=3
        ).with_read(
            "..................T.......C... ...C............................", sample_name=sn, n_fwd=3, n_rev=3
        ).with_output_phased_genotypes(True)

        svc_driver.call().with_output_vcf()\
            .has_record_for_variant(Variant("1", 18, "GTCAATGCATCCCATTG", "TTCAATGCCTCCCATTC"))\
            .with_sample(sn)\
            .has_genotype("1|1")
