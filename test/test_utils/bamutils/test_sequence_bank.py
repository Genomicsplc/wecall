# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from wecall.bamutils.sequence_bank import SequenceBank, AsciiVariantGenerator
from wecall.common.exceptions import EchidnaException
from wecall.genomics.reference_chromosome import ReferenceChromosome
from wecall.genomics.variant import Variant


class TestSequenceBank(TestCase):

    def test_should_fail_at_seq_with_different_length_to_reference(self):
        # Given
        ref_seq = "AAAA"
        seq = "CC"
        sequence_bank = SequenceBank(ReferenceChromosome(ref_seq))
        # Then
        self.assertRaises(EchidnaException, sequence_bank.add_sequence, seq)

    def test_should_be_able_to_add_snp_using_dsl_syntax(self):
        # Given
        input_ref = "CCC"
        snp_input = ".T."
        # When
        sequence_bank = SequenceBank(ReferenceChromosome(input_ref))
        sequence_bank.add_sequence(snp_input)
        read_lists = [builder.build_reads(0, {}) for builder in sequence_bank]
        reads = [read for read_list in read_lists for read in read_list]
        # Then
        self.assertEqual(reads[0].pos, 0)
        self.assertEqual(reads[0].seq, 'CTC')

    def test_should_be_able_to_add_snp_using_whitespace_dsl_syntax(self):
        # Given
        input_ref = "CC*AAGG"
        snp_input = "   .T. "
        # When
        sequence_bank = SequenceBank(ReferenceChromosome(input_ref))
        sequence_bank.add_sequence(snp_input)
        read_lists = [builder.build_reads(0, {}) for builder in sequence_bank]
        reads = [read for read_list in read_lists for read in read_list]
        # Then
        self.assertEqual(reads[0].pos, 2)
        self.assertEqual(reads[0].seq, 'ATG')


class TestAsciiVariantGenerator(TestCase):
    def test_should_generate_variant_from_ascii_text(self):
        ref = "ATAAAAAAAAAT"
        alt_1 = ".A........*."
        alt_2 = ".C.........."
        variant_generator = AsciiVariantGenerator(ReferenceChromosome(ref))

        gen_vars = variant_generator.get_variants([alt_1, alt_2])

        self.assertEqual(
            gen_vars,
            {
                Variant(variant_generator.reference.chrom, 1, "T", "A"),
                Variant(variant_generator.reference.chrom, 1, "T", "C"),
                Variant(variant_generator.reference.chrom, 9, "AA", "A")
            }
        )
