# All content Copyright (C) 2018 Genomics plc
import unittest
from wecall.genomics.chromosome import chromosome_comp


class TestChromosomeSort(unittest.TestCase):
    def test_should_sort_chrom_2_before_chrom_10(self):
        self.assertTrue(chromosome_comp('2', '10'))
        self.assertFalse(chromosome_comp('10', '2'))

    def test_should_sort_22_before_X(self):
        self.assertTrue(chromosome_comp('22', 'X'))
        self.assertFalse(chromosome_comp('X', '22'))

    def test_should_sort_X_before_Y(self):
        self.assertTrue(chromosome_comp('X', 'Y'))
        self.assertFalse(chromosome_comp('Y', 'X'))

    def test_should_sort_Y_before_MT(self):
        self.assertTrue(chromosome_comp('Y', 'MT'))
        self.assertFalse(chromosome_comp('MT', 'Y'))

    def test_should_sort_non_standard_chroms_after_standard_chroms(self):
        self.assertTrue(chromosome_comp("MT", "GL000193.1"))
        self.assertFalse(chromosome_comp("GL000193.1", "MT"))

    def test_should_sort_non_standard_chroms_lexicographically(self):
        self.assertTrue(chromosome_comp("GL000193.1", "GL000193.2"))
        self.assertFalse(chromosome_comp("GL000193.2", "GL000193.1"))
