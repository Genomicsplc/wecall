# All content Copyright (C) 2018 Genomics plc
import unittest

from wecall.genomics import variant
from wecall.utils.interval import ChromInterval


class TestvVariant(unittest.TestCase):

    def test_eq(self):
        reference = variant.Variant("1", 10, "ACGT", "ACC")
        self.assertEqual(reference, variant.Variant("1", 10, "ACGT", "ACC"))
        self.assertNotEqual(reference, variant.Variant("2", 10, "ACGT", "ACC"))
        self.assertNotEqual(reference, variant.Variant("1", 11, "ACGT", "ACC"))
        self.assertNotEqual(reference, variant.Variant("1", 10, "ACGC", "ACC"))
        self.assertNotEqual(reference, variant.Variant("1", 10, "ACGT", "ACT"))

    def test_ne(self):
        """
        Note: self.assert[Not]Equal() doesn't work in this case.
        """
        reference = variant.Variant("1", 10, "ACGT", "ACC")
        self.assertFalse(reference != variant.Variant("1", 10, "ACGT", "ACC"))
        self.assertTrue(reference != variant.Variant("2", 10, "ACGT", "ACC"))
        self.assertTrue(reference != variant.Variant("1", 11, "ACGT", "ACC"))
        self.assertTrue(reference != variant.Variant("1", 10, "ACGC", "ACC"))
        self.assertTrue(reference != variant.Variant("1", 10, "ACGT", "ACT"))

    def test_lt(self):
        reference = variant.Variant("1", 10, "ACGT", "ACC")
        self.assertFalse(reference < reference)
        self.assertGreaterEqual(
            reference, variant.Variant(
                "1", 10, "ACGT", "ACC"))
        self.assertLess(reference, variant.Variant("2", 10, "ACGT", "ACC"))
        self.assertLess(reference, variant.Variant("1", 11, "ACGT", "ACC"))
        self.assertLess(reference, variant.Variant("1", 10, "CCGT", "ACC"))
        self.assertLess(reference, variant.Variant("1", 10, "ACGT", "ACT"))

    def test_as_key(self):
        variant_obj = variant.Variant("1", 10, "ACGT", "ACC")
        expected_key = ("1", 10, 14, "ACGT", "ACC")
        self.assertEqual(
            hash(expected_key),
            hash(variant_obj),
            "Variant hash is broken.")

    def test_variant_type(self):
        chrom = "1"
        pos_from = 100
        self.assertEqual(variant.Variant(chrom, pos_from, 'T', '.').type, variant.TYPE_REF)
        self.assertEqual(variant.Variant(chrom, pos_from, 'C', 'C').type, variant.TYPE_REF)
        self.assertEqual(variant.Variant(chrom, pos_from, 'CTC', 'CTC').type, variant.TYPE_REF)
        self.assertEqual(variant.Variant(chrom, pos_from, 'T', 'TA').type, variant.TYPE_INS)
        self.assertEqual(variant.Variant(chrom, pos_from, 'CTG', 'CTAG').type, variant.TYPE_INS)
        self.assertEqual(variant.Variant(chrom, pos_from, 'TA', 'A').type, variant.TYPE_DEL)
        self.assertEqual(variant.Variant(chrom, pos_from, 'AT', 'A').type, variant.TYPE_DEL)
        self.assertEqual(variant.Variant(chrom, pos_from, 'T', 'A').type, variant.TYPE_SNP)
        self.assertEqual(variant.Variant(chrom, pos_from, 'TGT', 'AGT').type, variant.TYPE_SNP)
        self.assertEqual(variant.Variant(chrom, pos_from, 'TGTT', 'TGAT').type, variant.TYPE_SNP)
        self.assertEqual(variant.Variant(chrom, pos_from, 'AGTT', 'TGAT').type, variant.TYPE_MNP)
        self.assertEqual(variant.Variant(chrom, pos_from, 'AGTTATAT', 'TGATAAAT').type, variant.TYPE_MNP)
        self.assertEqual(variant.Variant(chrom, pos_from, 'A', '<INS:ME_ALU>').type, variant.TYPE_SYM)


class TestVariantOverlapsRegion(unittest.TestCase):

    def test_should_overlap_insertion_at_region_end(self):
        ins = variant.Variant('1', 15, 'A', 'AAAAA')
        self.assertTrue(ins.overlap(ChromInterval('1', 15, 16)))

    def test_should_not_overlap_insertion_after_region(self):
        ins = variant.Variant('1', 15, 'A', 'AAAAA')
        self.assertFalse(ins.overlap(ChromInterval('1', 14, 15)))

    def test_should_not_overlap_insertion_before_region(self):
        ins = variant.Variant('1', 15, 'A', 'AAAAA')
        self.assertFalse(ins.overlap(ChromInterval('1', 16, 20)))

    def test_should_overlap_deletion_at_region_start(self):
        deletion = variant.Variant('1', 15, 'AAAAA', 'A')
        self.assertTrue(deletion.overlap(ChromInterval('1', 19, 20)))

    def test_should_overlap_deletion_at_region_end(self):
        deletion = variant.Variant('1', 15, 'AAAAA', 'A')
        self.assertTrue(deletion.overlap(ChromInterval('1', 15, 16)))

    def test_should_not_overlap_deletion_before_region(self):
        deletion = variant.Variant('1', 15, 'AAAAA', 'A')
        self.assertFalse(deletion.overlap(ChromInterval('1', 14, 15)))

    def test_should_not_overlap_deletion_after_region(self):
        deletion = variant.Variant('1', 15, 'AAAAA', 'A')
        self.assertFalse(deletion.overlap(ChromInterval('1', 20, 21)))

    def test_should_overlap_mnp_at_region_start(self):
        mnp = variant.Variant('1', 15, 'ACGT', 'GCGC')
        self.assertTrue(mnp.overlap(ChromInterval('1', 18, 19)))

    def test_should_overlap_mnp_at_region_end(self):
        mnp = variant.Variant('1', 15, 'ACGT', 'GCGC')
        self.assertTrue(mnp.overlap(ChromInterval('1', 15, 16)))

    def test_should_not_overlap_mnp_before_region(self):
        mnp = variant.Variant('1', 15, 'ACGT', 'GCGC')
        self.assertFalse(mnp.overlap(ChromInterval('1', 14, 15)))

    def test_should_not_overlap_mnp_after_region(self):
        mnp = variant.Variant('1', 15, 'ACGT', 'GCGC')
        self.assertFalse(mnp.overlap(ChromInterval('1', 19, 20)))
