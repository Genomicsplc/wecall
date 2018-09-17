# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from wecall.bamutils.cigar import Cigar


class TestCigar(TestCase):
    def test_should_be_able_to_create_cigar_for_insertion(self):
        cigar = Cigar([(Cigar.INSERTION, 10)])
        self.assertEqual(str(cigar), "10I")

    def test_should_be_able_to_create_cigar_for_deletion(self):
        cigar = Cigar([(Cigar.DELETION, 10)])
        self.assertEqual(str(cigar), "10D")

    def test_should_be_able_to_create_cigar_for_match(self):
        cigar = Cigar([(Cigar.MATCH, 10)])
        self.assertEqual(str(cigar), "10M")

    def test_should_reduce_cigars_correctly_on_construction(self):
        cigar = Cigar([(Cigar.MATCH, 10), (Cigar.MATCH, 7)])
        self.assertEqual(str(cigar), "17M")

    def test_should_be_able_to_add_cigars(self):
        cigar_1 = Cigar([(Cigar.DELETION, 10), (Cigar.MATCH, 7)])
        cigar_2 = Cigar([(Cigar.MATCH, 6), (Cigar.INSERTION, 7)])

        self.assertEqual(str(cigar_1 + cigar_2), "10D13M7I")
