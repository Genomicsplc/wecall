# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.ascii_wecall_runner import AsciiWecallRunnerTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestVariantCallingFromTwoEqualLengthFwdReads(AsciiWecallRunnerTest):

    def test_calls_snp_with_one_one_genotype_on_forward_reads(self):
        self.calls_variants_with_genotype(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",  # input
            ["       ..............C.............       ",
             "  ...................C........            ",
             "    .................C..............      ",
             "    .................C..................  ",
             ".....................C......              "],

            [".....................C....................",  # expected output
             ".....................C...................."]
        )

    def test_calls_snp_with_zero_one_genotype_on_forward_reads(self):
        self.calls_variants_with_genotype(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ...................C........            ",
             "    ................................      ",
             "    .................C..................  ",
             ".....................C......              "],

            [".....................C....................",  # expected output
             ".........................................."]
        )

    def test_calls_two_snps_with_zero_one_genotypes(self):
        self.calls_variants_with_genotype(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       ..............C.............       ",
             "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,            ",
             "    ................................      ",
             "    ,,,,,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ",
             ".........G..................              "],

            [".....................C....................",  # expected output
             ".........G................................"]
        )

    def test_calls_snp_with_noise(self):
        self.calls_variants_with_genotype(
            "AAAAAAAAAAACGCACCCCCCATAAAAAAAATTTTTTTTTTT",
            ["       .....T........C.............       ",
             "            2                             ",
             "  ,,,,,,,,,,,,,,,,,,,c,,,,,,,,            ",
             "    ................................      ",
             "    ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,  ",
             "............................              "],

            [".....................C....................",  # expected output
             ".........................................."]
        )

    def test_calls_deletion_and_snp_at_same_location_in_repeat_region_with_few_reads_as_anchors(self):
        chrom = "1"
        sample = "sample"

        svc_driver = SVCDriver(self)
        svc_driver.with_ref_sequence(
            'CGAGAGAGAGAGAGAGAGAGATAGAGAGAGAGAGAGAGAGTC', chrom=chrom
        ).with_read(
            '....................**....................', n_rev=5, n_fwd=0, chrom=chrom, sample_name=sample
        ).with_read(
            '.....................G....................', n_rev=5, n_fwd=0, chrom=chrom, sample_name=sample
        )
        expect = svc_driver.call()
        vcf_expect = expect.with_output_vcf()
        vcf_expect \
            .has_record_for_variants(
                Variant(chrom, 21, "T", "G"),
                Variant(chrom, 19, "GAT", "G")
            ).with_sample(sample).has_phased_genotypes(".|1", "1|.")

    def test_calls_snp_as_hom_alt_call(self):
        # Based on real example found in NA12878 at SNP(2:15650796-15650797 G
        # --> A)
        self.calls_variants_with_genotype(
            #   7                                                                                                   8                                                                                                   9  # noqa
            # 9  0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0  # noqa
            # 789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567  # noqa
            "CATTTGCTCATTTGATTGCTTGAATTATCAAAATAGAACCTAGCCTATGCTGACTACTACAGCCCCAAGGACTGCCTGACACAAAATAGAAGCTCAATAAGTAAGTATAGAATGAATGGATGGATGGATGAATAGTCCATGTAAATATGGAACAACATGGAAATCTGGACTCTAATTCTGCCATTTGCTAGCTGGATGATC",  # noqa
            ["                                                                                                 ...A....................T.C......TT.T.CAAGCAG..G.CG.C.TACGAGATCGTGATGT....GG.G...A.A.G.G....CTTCC...C   ",  # noqa
             "                                                                                                 ...A.................................................................................................   ",  # noqa
             "                                                                                                 ...A.................................................................................................   ",  # noqa
             "                                                                                                 ...A.................................................................................................   ",  # noqa
             "                                                                                                 ...A.................................................................................................   ",  # noqa
             "                                                                                                 ...A.................................................................................................   ",  # noqa
             "                                                                                                 ...A.................................................................................................   ",  # noqa
             "                                                                                                 ...A.................................................................................................   ", ],  # noqa

            ["....................................................................................................A....................................................................................................",  # noqa
             "....................................................................................................A...................................................................................................."]  # noqa
        )
