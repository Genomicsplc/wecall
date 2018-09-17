# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestCallingForTrio(AsciiWecallRunnerTest):

    def test_multi_sample_variant_calling(self):
        reference = "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT"
        samples = {
            "NA12878": ["       ..............C.............       ",
                        "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,           ",
                        "    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,     ",
                        "..............*.............              "],

            "NA12891": ["       ..............C.............       ",
                        "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,           ",
                        "    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,     ",
                        "..............*.............              "],

            "NA12892": ["       ..............C.............       ",
                        "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,           ",
                        "    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,     ",
                        "..............*.............              "]
        }

        expected_haplotypes = {
            "NA12878": ["..............*...........................",
                        ".....................C...................."],

            "NA12891": ["..............*...........................",
                        ".....................C...................."],

            "NA12892": ["..............*...........................",
                        ".....................C...................."]
        }

        self.calls_variants_from_samples(
            reference, samples, expected_haplotypes)

    def test_multi_sample_variant_calling_with_some_regions_covered_in_only_one_sample_new(self):
        chrom = "1"
        driver = SVCDriver(self)
        driver.with_ref_sequence(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTTCTCTAATTGGGTCACGTATGCATGACGTTGTGGGGAACCCCTGG", chrom=chrom
        ).with_read(
            "       ..............C.............                                                   ", sample_name="NA12878"  # noqa
        ).with_read(
            "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,                                                       ", sample_name="NA12878"  # noqa
        ).with_read(
            "    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,                                                 ", sample_name="NA12878"  # noqa
        ).with_read(
            "..............*.............                                                          ", sample_name="NA12878"  # noqa
        ).with_read(  # noqa
            "                                                       ..............A................", sample_name="NA12878"  # noqa
        ).with_read(
            "                                                       ..............A................", sample_name="NA12878"  # noqa
        ).with_read(
            "                                                       ..............A................", sample_name="NA12878"  # noqa
        ).with_read(
            "                                                       ..............A................", sample_name="NA12878"  # noqa
        ).with_read(
            "       ..............C.............                                                   ", sample_name="NA12891"  # noqa
        ).with_read(
            "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,                                                       ", sample_name="NA12891"  # noqa
        ).with_read(
            "    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,                                                 ", sample_name="NA12891"  # noqa
        ).with_read(
            "..............*.............                                                          ", sample_name="NA12891"  # noqa
        ).with_read(
            "       ..............C.............                                                   ", sample_name="NA12892"  # noqa
        ).with_read(
            "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,,                                                       ", sample_name="NA12892"  # noqa
        ).with_read(
            "    ,,,,,,,,,,*,,,,,,,,,,,,,,,,,,,,,,                                                 ", sample_name="NA12892"  # noqa
        ).with_read(
            "..............*.............                                                          ", sample_name="NA12892"  # noqa
        )

        expect = driver.call()

        vcf = expect \
            .with_output_vcf() \
            .record_count(3)

        C_CA = vcf.has_record_for_variant(Variant(chrom, 13, "CA", "C"))
        C_CA.with_sample("NA12878").has_genotype("0/1")
        C_CA.with_sample("NA12891").has_genotype("0/1")
        C_CA.with_sample("NA12892").has_genotype("0/1")

        A_C = vcf.has_record_for_variant(Variant(chrom, 21, "A", "C"))
        A_C.with_sample("NA12878").has_genotype("0/1")
        A_C.with_sample("NA12891").has_genotype("0/1")
        A_C.with_sample("NA12892").has_genotype("0/1")

        T_A = vcf.has_record_for_variant(Variant(chrom, 69, "T", "A"))
        T_A.with_sample("NA12878").has_genotype("1/1")
        T_A.with_sample("NA12891").has_genotype("./.")
        T_A.with_sample("NA12892").has_genotype("./.")
