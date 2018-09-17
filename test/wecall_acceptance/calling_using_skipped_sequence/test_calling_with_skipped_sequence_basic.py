# All content Copyright (C) 2018 Genomics plc
from wecall.genomics.variant import Variant
from wecall_test_drivers.base_test import BaseTest
from wecall_test_drivers.svc_driver import SVCDriver


class TestCallingUsingSkippedSequenceBasic(BaseTest):
    def setUp(self):
        BaseTest.setUp(self)
        self.scv = SVCDriver(self)
        self.scv.with_turn_on_large_variant_calls(True)
        self.scv.with_verbosity(6)

    def test_should_call_variants_minimal_example(self):
        self.scv.with_verbosity(6)
        self.scv.with_ref_sequence(
            "ATAAAAAATATGTACATAAAAATCAAAATCAAAGAAAGAACATGCAGTAGCTGAAAAAAAATATCTTCTCACCCTAAAACTGCTCTATGTTTTAAACTATTATTGCTAGGATCACTAGGACTTAGTAAAAAGCAATGCCTTACACAGGCAAC",  # noqa
            pos_from=10
        ).with_read(
            "                                                        ..............***********........                                                               ",  # noqa
            cigar_string="14M8S", read_mate_start=91, n_rev=2, n_fwd=2, read_start=66
        ).with_read(
            "                                                               .......***********..............                                                         ",  # noqa
            cigar_string="7S14M", read_start=91, n_rev=2, n_fwd=2, read_mate_start=66
        )
        expect = self.scv.call()
        expect.with_output_vcf() \
            .has_record_for_variant(Variant("1", 79, 'CACCCTAAAACT', 'C'))

    def test_should_call_variants_simple_example(self):
        l_ref_padding = "CTTAAAGTGTAATAAAAAATATGTACATAAAAATCAAAATCAAAGAAAGAACATGCAGTAGCTGAAAAAAAATATCTTCTC"
        r_ref_padding = "GCTCTATGTTTTAAACTATTATTGCTAGGATCACTAGGACTTAGTAAAAAGCAATGCCTTACACAGGCAACAA"
        overlap = "AAAAGCAT"
        remaining_deleted = "ACCCTAAAACTTAGAGTATTCTCAATAAAAAAAAAAAAATTAAAAAAAAA"
        seq_buffer = 4  # amount of soft clipping
        self.call_deletion(l_ref_padding, overlap, remaining_deleted, r_ref_padding, seq_buffer)

    def test_should_call_variants_simple_example_no_overlap(self):
        l_ref_padding = "CTTAAAGTGTAATAAAAAATATGTACATAAAAATCAAAATCAAAGAAAGAACATGCAGTAGCTGAAAAAAAATATCTTCTC"
        r_ref_padding = "GCTCTATGTTTTAAACTATTATTGCTAGGATCACTAGGACTTAGTAAAAAGCAATGCCTTACACAGGCAACAA"
        overlap = ""
        remaining_deleted = "ACCCTAAAACTTAGAGTATTCTCAATAAAAAAAAAAAAATTAAAAAAAAA"
        seq_buffer = 8  # amount of soft clipping
        self.call_deletion(l_ref_padding, overlap, remaining_deleted, r_ref_padding, seq_buffer)

    def test_should_call_variants_simple_example_some_overlap(self):
        l_ref_padding = "CTTAAAGTGTAATAAAAAATATGTACATAAAAATCAAAATCAAAGAAAGAACATGCAGTAGCTGAAAAAAAATATCTTCTC"
        r_ref_padding = "GCTCTATGTTTTAAACTATTATTGCTAGGATCACTAGGACTTAGTAAAAAGCAATGCCTTACACAGGCAACAA"
        overlap = "GCAT"
        remaining_deleted = "ACCCTAAAACTTAGAGTATTCTCAATAAAAAAAAAAAAATTAAAAAAAAA"
        seq_buffer = 6  # amount of soft clipping
        self.call_deletion(l_ref_padding, overlap, remaining_deleted, r_ref_padding, seq_buffer)

    def call_deletion(self, l_ref_padding, overlap, remaining_deleted, r_ref_padding, seq_buffer):
        chrom = '20'
        ref_start = 9696680
        sample_name = 'sample'

        reference_sequence = l_ref_padding + overlap + remaining_deleted + overlap + r_ref_padding

        del_str = "*" * (len(overlap) + len(remaining_deleted))

        alt_read_1 = l_ref_padding + overlap + del_str + "." * seq_buffer + " " * (len(r_ref_padding) - seq_buffer)
        ref_read_1 = l_ref_padding + overlap + "." * seq_buffer + " " * (len(r_ref_padding) - seq_buffer + len(del_str))

        alt_read_2 = " " * (len(l_ref_padding) - seq_buffer) + "." * seq_buffer + del_str + overlap + r_ref_padding
        ref_read_2 = " " * (len(l_ref_padding) - seq_buffer + len(del_str)) + "." * seq_buffer + overlap + r_ref_padding

        event_start = ref_start + len(l_ref_padding)
        event_end = event_start + len(overlap) + len(remaining_deleted)

        self.scv.with_ref_sequence(
            reference_sequence, chrom=chrom, pos_from=ref_start
        ).with_read(
            alt_read_1, n_fwd=2, n_rev=2, chrom=chrom, sample_name=sample_name, read_start=ref_start,
            cigar_string='{}M{}S'.format(len(l_ref_padding) + len(overlap), seq_buffer), read_mate_start=event_end
        ).with_read(
            alt_read_2, n_fwd=2, n_rev=2, chrom=chrom, sample_name=sample_name, read_start=event_end,
            cigar_string='{}S{}M'.format(seq_buffer, len(overlap) + len(r_ref_padding)), read_mate_start=ref_start
        ).with_read(
            ref_read_1, n_fwd=2, n_rev=2, chrom=chrom, sample_name=sample_name, read_start=ref_start,
            cigar_string='{}M'.format(len(l_ref_padding) + len(overlap) + seq_buffer), read_mate_start=event_end
        ).with_read(
            ref_read_2, n_fwd=2, n_rev=2, chrom=chrom, sample_name=sample_name, read_start=event_end - seq_buffer,
            cigar_string='{}M'.format(seq_buffer + len(overlap) + len(r_ref_padding)), read_mate_start=ref_start)

        expect = self.scv.call()

        base_before = reference_sequence[len(l_ref_padding) - 1:len(l_ref_padding)]

        variant = Variant(
            chrom,
            event_start - 1,
            base_before + overlap + remaining_deleted,
            base_before
        )

        expect.with_output_vcf()\
            .has_record_for_variant(variant)\
            .with_sample(sample_name)\
            .has_genotype('0/1')
