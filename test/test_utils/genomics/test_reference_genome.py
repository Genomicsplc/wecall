# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from wecall.genomics.reference_genome import InMemoryReferenceGenome


class TestInMemoryReferenceGenome(TestCase):

    def test_should_not_get_duplicate_or_extra_chromosomes(self):
        reference_genome = InMemoryReferenceGenome()
        reference_genome.with_chrom("1", "", 0)
        reference_genome.with_chrom("1", "", 0)
        reference_genome.with_chrom("2", "", 0)

        self.assertEqual({"1", "2"}, set(reference_genome.chromosomes()))

    def test_should_get_chromosomes_in_correct_order(self):
        reference_genome = InMemoryReferenceGenome()
        reference_genome.with_chrom("1", "", 0)
        reference_genome.with_chrom("3", "", 0)
        reference_genome.with_chrom("2", "", 0)

        self.assertEqual(['1', '2', '3'], reference_genome.chromosomes())

    def test_should_get_total_chromosome_length(self):
        reference_genome = InMemoryReferenceGenome()
        reference_genome.with_chrom("1", "ATG", 11)

        self.assertEqual(14, reference_genome.get_chrom_length("1"))

    def test_should_fetch_correct_sequence_with_padding(self):
        reference_genome = InMemoryReferenceGenome()
        reference_genome.with_chrom("1", "ATG", 11)

        self.assertEqual("NA", reference_genome.fetch("1", 10, 12))

    def test_should_get_index_error_when_out_of_range(self):
        reference_genome = InMemoryReferenceGenome()
        reference_genome.with_chrom("1", "ATG", 11)

        with self.assertRaises(IndexError):
            print((reference_genome.fetch("1", 100, 120)))
