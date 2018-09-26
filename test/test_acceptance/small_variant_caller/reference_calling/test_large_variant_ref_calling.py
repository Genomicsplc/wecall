from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestCallingLargeVariantWithReferenceCalling(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.scv = SVCDriver(self)
        self.scv.with_turn_on_large_variant_calls(True)
        self.scv.with_output_ref_calls(True)
        self.scv.with_verbosity(6)

    def test_should_call_reference_for_minimal_large_variant_example(self):
        self.scv.with_ref_sequence(
            "ATAAAAAATATGTACATAAAAATCAAAATCAAAGAAAGAACATGCAGTAGCTGAAAAAAAATATCTTCTCACCCTAAAACTGCTCTATGTTTTAAACTATTATTGCTAGGATCACTAGGACTTAGTAAAAAGCAATGCCTTACACAGGCAAC",  # noqa
            pos_from=10
        ).with_read(
            "                                                        ..............***********........                                                               ",  # noqa
            cigar_string="14M8S", read_mate_start=91, n_rev=10, n_fwd=10, read_start=66
        ).with_read(
            "                                                               .......***********..............                                                         ",  # noqa
            cigar_string="7S14M", read_start=91, n_rev=10, n_fwd=10, read_mate_start=66
        )

        expect = self.scv.call().with_output_vcf()

        expect.has_record_for_variant(Variant("1", 79, 'CACCCTAAAACT', 'C'))
        expect.has_reference_calls_for_region("1", 60, 66)  # first read starting
        expect.has_reference_calls_for_region("1", 66, 79)  # second read starting
        expect.has_reference_calls_for_region("1", 91, 105)  # after large deletion

    def test_should_call_reference_for_minimal_large_variant_example_for_two_samples(self):
        self.scv.with_ref_sequence(
            "ATAAAAAATATGTACATAAAAATCAAAATCAAAGAAAGAACATGCAGTAGCTGAAAAAAAATATCTTCTCACCCTAAAACTGCTCTATGTTTTAAACTATTATTGCTAGGATCACTAGGACTTAGTAAAAAGCAATGCCTTACACAGGCAAC",  # noqa
            pos_from=10
        ).with_read(
            "                                                        ..............***********........                                                               ",  # noqa
            cigar_string="14M8S", read_mate_start=91, n_rev=10, n_fwd=10, read_start=66, sample_name='sample1'
        ).with_read(
            "                                                               .......***********..............                                                         ",  # noqa
            cigar_string="7S14M", read_start=91, n_rev=10, n_fwd=10, read_mate_start=66, sample_name='sample2'
        ).with_read(
            "                                                        ...............................                                                                 ",  # noqa
            cigar_string="31M", n_rev=10, n_fwd=10, read_start=66, sample_name='sample1'
        ).with_read(
            "                                                               ......................G.......                                                           ",  # noqa
            cigar_string="30M", read_start=73, n_rev=20, n_fwd=20, sample_name='sample1'
        )

        expect = self.scv.call().with_output_vcf()

        expect.has_record_for_variant(Variant("1", 79, 'CACCCTAAAACT', 'C')).with_sample('sample1').has_genotype('1|0')
        expect.has_record_for_variant(Variant("1", 95, 'T', 'G')).with_sample('sample1').has_genotype('0|1')

        expect.has_record_for_variant(Variant("1", 79, 'CACCCTAAAACT', 'C')).with_sample('sample2').has_genotype('1|1')
        expect.has_record_for_variant(Variant("1", 95, 'T', 'G')).with_sample('sample2').has_genotype('0|0')

        # only get reference calls when there is read data
        expect.has_reference_calls_for_region("1", 73, 79)
        expect.has_reference_calls_for_region("1", 91, 95)
        expect.has_reference_calls_for_region("1", 96, 97)
