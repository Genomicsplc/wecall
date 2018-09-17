# All content Copyright (C) 2018 Genomics plc
from unittest import TestCase
from pysam import AlignedRead
from wecall.bamutils.bam_region_iterator import BAMRegionIterator
from wecall.utils.interval import Interval


def make_read(pos, aend):
    read = AlignedRead()
    read.pos = pos
    read.seq = 'N' * (aend - pos)
    read.cigarstring = '{}M'.format(aend - pos)
    assert read.aend == aend, '{} != {}'.format(read.aend, aend)
    return read


class TestBAMFileIterator(TestCase):
    def test_should_return_all_reads_if_entire_fetch_region_requested(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        actual_results = list(iterator(Interval(0, 4)))
        self.assertEqual(all_reads, actual_results)

    def test_should_only_return_reads_with_start_overlapping_interval(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        actual_results = list(iterator(Interval(0, 2)))
        self.assertEqual(all_reads[0:2], actual_results)

    def test_should_only_return_reads_with_end_overlapping_interval(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        actual_results = list(iterator(Interval(2, 4)))
        self.assertEqual(all_reads[2:4], actual_results)

    def test_should_return_overlapping_reads_from_middle_of_query_region(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
            make_read(4, 5),
            make_read(5, 6),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        actual_results = list(iterator(Interval(2, 4)))
        self.assertEqual(all_reads[2:4], actual_results)

    def test_should_return_nothing_if_interval_does_not_overlap_any_reads(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        actual_results = list(iterator(Interval(4, 6)))
        self.assertEqual([], actual_results)

    def test_should_return_two_sets_of_values_for_intervals_in_correct_order(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
            make_read(4, 5),
            make_read(5, 6),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        self.assertEqual(all_reads[0:2], list(iterator(Interval(0, 2))))
        self.assertEqual(all_reads[3:5], list(iterator(Interval(3, 5))))

    def test_should_not_miss_any_reads_between_regions(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
            make_read(4, 5),
            make_read(5, 6),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        self.assertEqual(all_reads[0:2], list(iterator(Interval(0, 2))))
        self.assertEqual(all_reads[2:3], list(iterator(Interval(2, 3))))

    def test_should_fail_for_regions_in_incorrect_order(self):
        all_reads = [
            make_read(0, 1),
            make_read(1, 2),
            make_read(2, 3),
            make_read(3, 4),
            make_read(4, 5),
            make_read(5, 6),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        self.assertEqual(all_reads[3:5], list(iterator(Interval(3, 5))))
        with self.assertRaises(Exception):
            print((iterator(Interval(0, 2))))

    def test_should_return_reads_that_overlap_both_intervals_twice(self):
        all_reads = [
            make_read(0, 2),
            make_read(1, 3),
        ]

        iterator = BAMRegionIterator((read for read in all_reads))

        self.assertEqual(all_reads[0:1], list(iterator(Interval(0, 1))))
        self.assertEqual(all_reads[0:2], list(iterator(Interval(1, 2))))
