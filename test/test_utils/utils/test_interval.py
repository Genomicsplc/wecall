# All content Copyright (C) 2018 Genomics plc
from wecall.common.exceptions import EchidnaException
from wecall.utils import interval
import unittest
from wecall.utils.interval import read_chrom_intervals, read_chrom_interval, read_interval, ChromInterval, \
    merge_intervals


class TestInterval(unittest.TestCase):

    def setUp(self):
        ref_start = 10
        ref_end = 40
        ref_mid = (ref_start + ref_end) // 2

        test_pos = [
            ref_start - 5,
            ref_start - 1,
            ref_start,
            ref_mid,
            ref_end - 1,
            ref_end,
            ref_end + 5,
        ]

        self.ref_interval = interval.Interval(ref_start, ref_end)

        # Matrix of expected output from contigIntervalsOverlap function. Deals
        # explicitly with corner cases instead of relying on alternative
        # function.
        overlaps = [
            [False, False, True, True, True, True, True, ],
            [False, False, True, True, True, True, True, ],
            [False, False, True, True, True, True, True, ],
            [False, False, False, True, True, True, True, ],
            [False, False, False, False, True, True, True, ],
            [False, False, False, False, False, True, True, ],
            [False, False, False, False, False, False, False, ],
        ]

        valid_intervals = [
            [True, True, True, True, True, True, True, ],
            [False, True, True, True, True, True, True, ],
            [False, False, True, True, True, True, True, ],
            [False, False, False, True, True, True, True, ],
            [False, False, False, False, True, True, True, ],
            [False, False, False, False, False, True, True, ],
            [False, False, False, False, False, False, True, ],
        ]

        intersects = [
            [None, None, (ref_start, ref_start), (ref_start, test_pos[3]), (ref_start, test_pos[4]), (ref_start, ref_end), (ref_start, ref_end), ],  # noqa
            [None, None, (ref_start, ref_start), (ref_start, test_pos[3]), (ref_start, test_pos[4]), (ref_start, ref_end), (ref_start, ref_end), ],  # noqa
            [None, None, (ref_start, ref_start), (ref_start, test_pos[3]), (ref_start, test_pos[4]), (ref_start, ref_end), (ref_start, ref_end), ],  # noqa
            [None, None, None, (test_pos[3], test_pos[3]), (test_pos[3], test_pos[4]), (test_pos[3], ref_end), (test_pos[3], ref_end), ],  # noqa
            [None, None, None, None, (test_pos[4], test_pos[4]), (test_pos[4], ref_end), (test_pos[4], ref_end), ],  # noqa
            [None, None, None, None, None, (ref_end, ref_end), (ref_end, ref_end), ],  # noqa
            [None, None, None, None, None, None, None, ],  # noqa
        ]

        combinations = [
            [(ref_start, ref_end), (ref_start, ref_end), (test_pos[0], ref_end), (test_pos[0], ref_end), (test_pos[0], ref_end), (test_pos[0], ref_end), (test_pos[0], test_pos[6]), ],  # noqa
            [(ref_start, ref_end), (ref_start, ref_end), (test_pos[1], ref_end), (test_pos[1], ref_end), (test_pos[1], ref_end), (test_pos[1], ref_end), (test_pos[1], test_pos[6]), ],  # noqa
            [(ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, test_pos[6]), ],  # noqa
            [(ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, test_pos[6]), ],  # noqa
            [(ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, test_pos[6]), ],  # noqa
            [(ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, test_pos[6]), ],  # noqa
            [(ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), ],  # noqa
        ]

        left_differences = [
            [(ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), ],  # noqa
            [None, (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), ],  # noqa
            [None, None, (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), (ref_start, ref_start), ],  # noqa
            [None, None, None, (ref_start, test_pos[3]), (ref_start, test_pos[3]), (ref_start, test_pos[3]), (ref_start, test_pos[3]), ],  # noqa
            [None, None, None, None, (ref_start, test_pos[4]), (ref_start, test_pos[4]), (ref_start, test_pos[4]), ],  # noqa
            [None, None, None, None, None, (ref_start, ref_end), (ref_start, ref_end), ],  # noqa
            [None, None, None, None, None, None, (ref_start, ref_end), ],  # noqa
        ]

        right_differences = [
            [(ref_start, ref_end), (ref_start, ref_end), (ref_start, ref_end), (test_pos[3], ref_end), (test_pos[4], ref_end), (ref_end, ref_end), (ref_end, ref_end), ],  # noqa
            [None, (ref_start, ref_end), (ref_start, ref_end), (test_pos[3], ref_end), (test_pos[4], ref_end), (ref_end, ref_end), (ref_end, ref_end), ],  # noqa
            [None, None, (ref_start, ref_end), (test_pos[3], ref_end), (test_pos[4], ref_end), (ref_end, ref_end), (ref_end, ref_end), ],  # noqa
            [None, None, None, (test_pos[3], ref_end), (test_pos[4], ref_end), (ref_end, ref_end), (ref_end, ref_end), ],  # noqa
            [None, None, None, None, (test_pos[4], ref_end), (ref_end, ref_end), (ref_end, ref_end), ],  # noqa
            [None, None, None, None, None, (ref_end, ref_end), (ref_end, ref_end), ],  # noqa
            [None, None, None, None, None, None, (ref_end, ref_end), ],  # noqa
        ]

        self.test_intervals = []

        # Expecting no exceptions.
        for start_index, start_pos in enumerate(test_pos):
            for end_index, end_pos in enumerate(test_pos):

                test_interval = interval.Interval(start_pos, end_pos)

                intersect = self.ref_interval.intersection(test_interval)
                combination = self.ref_interval.combination(test_interval)

                left_difference = self.ref_interval.leftDifference(test_interval)  # noqa
                right_difference = self.ref_interval.rightDifference(test_interval)  # noqa
                union = left_difference.combination(intersect).combination(right_difference)  # noqa

                # Find actual overlap with interval.
                actual_overlap = test_interval.overlapOrTouch(self.ref_interval)  # noqa
                # Find actual overlap, reverse arguments.
                ractual_overlap = self.ref_interval.overlapOrTouch(test_interval)  # noqa

                # Note in tests use expected_overlap rather than actual_overlap or ractual_overlap.  # noqa
                self.test_intervals.append({
                    "interval": test_interval,
                    "expected_overlap": overlaps[start_index][end_index],
                    "actual_overlap": actual_overlap,
                    "ractual_overlap": ractual_overlap,
                    "start_index": start_index,
                    "end_index": end_index,
                    "expected_intersect": intersects[start_index][end_index],
                    "actual_intersect": intersect,
                    "expected_combination": combinations[start_index][end_index],  # noqa
                    "actual_combination": combination,
                    "expected_left_difference": left_differences[start_index][end_index],  # noqa
                    "actual_left_difference": left_difference,
                    "expected_right_difference": right_differences[start_index][end_index],  # noqa
                    "actual_right_difference": right_difference,
                    "actual_union": union,
                    "valid_interval": valid_intervals[start_index][end_index],
                })

        # Manually add interesting cases.

    def test_contigIntervalsOverlap(self):

        for item in self.test_intervals:
            self.assertTrue(
                item["expected_overlap"] is item["actual_overlap"],
                ", ".join((x.format(self.ref_interval, **item) for x in (
                    # test-specific data
                    "expected: {expected_overlap}",
                    "actual: {actual_overlap}",
                    # data common to all tests
                    "interval: {interval!r}",
                    "reference: {0!r}",
                    "start-index: {start_index}",
                    "end-index: {end_index}",
                )))
            )
            self.assertTrue(
                item["expected_overlap"] is item["ractual_overlap"],
                ", ".join((x.format(self.ref_interval, **item) for x in (
                    # test-specific data
                    "expected: {expected_overlap}",
                    "r-actual: {ractual_overlap}",
                    # data common to all tests
                    "interval: {interval!r}",
                    "reference: {0!r}",
                    "start-index: {start_index}",
                    "end-index: {end_index}",
                )))
            )

    def __standardTest(self, expected, actual):
        return actual.is_null \
            if expected is None else actual is not None and actual.start == expected[0] and actual.end == expected[1]

    def test_intersection(self):

        for item in self.test_intervals:
            expected = item["expected_intersect"]
            actual = item["actual_intersect"]

            self.assertTrue(
                self.__standardTest(expected, actual),
                ", ".join((x.format(self.ref_interval, **item) for x in (
                    # test-specific data
                    "expected: {expected_intersect}",
                    "actual: {actual_intersect}",
                    # data common to all tests
                    "interval: {interval!r}",
                    "reference: {0!r}",
                    "start-index: {start_index}",
                    "end-index: {end_index}",
                )))
            )

    def test_combine(self):
        for item in self.test_intervals:
            expected = item["expected_combination"]
            actual = item["actual_combination"]

            self.assertTrue(
                self.__standardTest(expected, actual),
                ", ".join((x.format(self.ref_interval, **item) for x in (
                    # test-specific data
                    "expected: {expected_combination}",
                    "actual: {actual_combination}",
                    # data common to all tests
                    "interval: {interval!r}",
                    "reference: {0!r}",
                    "start-index: {start_index}",
                    "end-index: {end_index}",
                )))
            )

    def test_left_differences(self):
        for item in self.test_intervals:
            expected = item["expected_left_difference"]
            actual = item["actual_left_difference"]

            self.assertTrue(
                self.__standardTest(expected, actual),
                ", ".join((x.format(self.ref_interval, **item) for x in (
                    # test-specific data
                    "expected: {expected_left_difference}",
                    "actual: {actual_left_difference}",
                    # data common to all tests
                    "interval: {interval!r}",
                    "reference: {0!r}",
                    "start-index: {start_index}",
                    "end-index: {end_index}",
                )))
            )

    def test_right_differences(self):
        for item in self.test_intervals:
            expected = item["expected_right_difference"]
            actual = item["actual_right_difference"]

            self.assertTrue(
                self.__standardTest(expected, actual),
                ", ".join((x.format(self.ref_interval, **item) for x in (
                    # test-specific data
                    "expected: {expected_right_difference}",
                    "actual: {actual_right_difference}",
                    # data common to all tests
                    "interval: {interval!r}",
                    "reference: {0!r}",
                    "start-index: {start_index}",
                    "end-index: {end_index}",
                )))
            )

    def test_partition_union_validity(self):
        """
        Test validity of the result of `L(A\B) + (AnB) + R(A\B)`.
        """

        for item in self.test_intervals:

            failure_message = ", ".join((x.format(self.ref_interval, **item) for x in (
                # test-specific data
                "left_difference: {actual_left_difference}",
                "intersect: {actual_intersect}",
                "right_difference: {actual_right_difference}",
                "union: {actual_union}",
                # data common to all tests
                "interval: {interval!r}",
                "reference: {0!r}",
                "start-index: {start_index}",
                "end-index: {end_index}",
            )))

            self.assertTrue(
                item["actual_union"].is_null != item["valid_interval"],
                failure_message
            )

    def test_partition_union(self):
        """
        Tests that `A = L(A\B) + (AnB) + R(A\B)` when B is valid.
        """

        for item in self.test_intervals:

            failure_message = ", ".join((x.format(self.ref_interval, **item) for x in (
                # test-specific data
                "left_difference: {actual_left_difference}",
                "intersect: {actual_intersect}",
                "right_difference: {actual_right_difference}",
                "union: {actual_union}",
                # data common to all tests
                "interval: {interval!r}",
                "reference: {0!r}",
                "start-index: {start_index}",
                "end-index: {end_index}",
            )))

            if item["valid_interval"]:

                self.assertTrue(
                    self.ref_interval == item["actual_union"],
                    failure_message
                )

    def test_contigInterval_equal(self):
        contigs = ["chr1", "chr2"]
        starts = [0, 5]
        ends = [10, 20]
        ref_interval = interval.Interval(starts[0], ends[0])

        self.assertEqual(ref_interval, ref_interval)

        for index1, contig in enumerate(contigs):
            for index2, start in enumerate(starts):
                for index3, end in enumerate(ends):
                    if index2 == 0 and index3 == 0:
                        self.assertEqual(
                            ref_interval,
                            interval.Interval(start, end)
                        )
                    else:
                        self.assertNotEqual(
                            ref_interval,
                            interval.Interval(start, end)
                        )

    def test_interval_serialisation(self):
        interval_a = interval.Interval(10, 20)
        interval_b = interval.Interval(30, 40)

        self.assertEqual(
            interval_a,
            interval.Interval.fromDict(interval.Interval.toDict(interval_a))
        )

        self.assertEqual(
            interval_b,
            interval.Interval.fromDict(interval.Interval.toDict(interval_b))
        )


class TestIntervalCoverage(unittest.TestCase):

    def test_overlap_size(self):
        self.assertEqual(
            interval.Interval.overlap_size(
                interval.Interval(0, 10),
                interval.Interval(20, 300),
            ),
            None
        )

        self.assertEqual(
            interval.Interval.overlap_size(
                interval.Interval(0, 10),
                interval.Interval(5, 12),
            ),
            5
        )

        self.assertEqual(
            interval.Interval.overlap_size(
                interval.Interval(0, 10),
                interval.Interval(2, 5),
            ),
            3
        )

    def test_interval_coverage_zero_length(self):
        self.assertEqual(
            interval.interval_coverage(
                {}, interval.Interval(
                    0, 10)), 0)

        self.assertEqual(
            interval.interval_coverage(
                {interval.Interval(1, 1)},
                interval.Interval(0, 10)
            ),
            0
        )

    def test_interval_coverage_outside_main_interval(self):
        self.assertEqual(
            interval.interval_coverage(
                {
                    interval.Interval(-2, 0),
                    interval.Interval(10, 12)
                },
                interval.Interval(0, 10)
            ),
            0
        )

    def test_interval_coverage_overlap_left(self):
        self.assertEqual(
            interval.interval_coverage(
                {interval.Interval(-2, 1)},
                interval.Interval(0, 10)
            ),
            1
        )

    def test_interval_coverage_overlap_right(self):
        self.assertEqual(
            interval.interval_coverage(
                {interval.Interval(9, 11)},
                interval.Interval(0, 10)
            ),
            1
        )

    def test_interval_coverage_inclusion(self):
        self.assertEqual(
            interval.interval_coverage(
                {interval.Interval(-1, 11)},
                interval.Interval(0, 10)
            ),
            10
        )

        self.assertEqual(
            interval.interval_coverage(
                {interval.Interval(0, 10)},
                interval.Interval(0, 10)
            ),
            10
        )

    def test_interval_coverage_middle(self):
        self.assertEqual(
            interval.interval_coverage(
                {
                    interval.Interval(1, 4),
                    interval.Interval(5, 6)
                },
                interval.Interval(0, 10)
            ),
            4
        )

        self.assertEqual(
            interval.interval_coverage(
                {
                    interval.Interval(5, 6),
                    interval.Interval(1, 4),
                },
                interval.Interval(0, 10)
            ),
            4
        )

        self.assertEqual(
            interval.interval_coverage(
                {
                    interval.Interval(1, 4),
                    interval.Interval(3, 6)
                },
                interval.Interval(0, 10)
            ),
            5
        )

        self.assertEqual(
            interval.interval_coverage(
                {
                    interval.Interval(1, 3),
                    interval.Interval(1, 6),
                    interval.Interval(5, 7),
                },
                interval.Interval(0, 10)
            ),
            6
        )

        self.assertEqual(
            interval.interval_coverage(
                {
                    interval.Interval(1, 6),
                    interval.Interval(1, 3),
                    interval.Interval(5, 7),
                },
                interval.Interval(0, 10)
            ),
            6
        )


class TestChromInterval(unittest.TestCase):
    def test_equality(self):
        self.assertTrue(ChromInterval("1", 0, 0) == ChromInterval("1", 0, 0))
        self.assertFalse(ChromInterval("1", 0, 0) != ChromInterval("1", 0, 0))

    def test_overlap(self):
        intvl1 = interval.ChromInterval("chr1", 1, 10)
        intvl2 = interval.ChromInterval("chr2", 1, 10)
        self.assertFalse(intvl1.overlap(intvl2))
        self.assertTrue(intvl1.overlap(intvl1))

    def test_should_raise_assertion_for_different_chroms(self):
        intvl1 = interval.ChromInterval("chr1", 1, 10)
        intvl2 = interval.ChromInterval("chr2", 1, 10)
        self.assertRaises(
            AssertionError,
            interval.ChromInterval.intersection,
            intvl1,
            intvl2)

    def test_should_form_sensible_chrom_interval(self):
        intvl1 = interval.ChromInterval("chr1", 1, 10)
        intvl2 = interval.ChromInterval("chr1", 5, 15)
        self.assertEqual(
            intvl1.intersection(intvl2),
            interval.ChromInterval(
                "chr1",
                5,
                10))

    def test_chrom_interval_in_short_form_should_overlap(self):
        intvl1 = interval.ChromInterval("chr1")
        intvl2 = interval.ChromInterval("chr1", 5, 15)
        self.assertTrue(intvl1.overlap(intvl2))

    def test_chrom_interval_in_short_form_should_give_valid_intersection(self):
        intvl1 = interval.ChromInterval("chr1")
        intvl2 = interval.ChromInterval("chr1", 5, 15)
        self.assertEqual(intvl1.intersection(intvl2), intvl2)

    def test_should_initialise_from_region_string(self):
        interval = ChromInterval.from_string("20:1-10")
        self.assertEqual(interval.chrom, "20")
        self.assertEqual(interval.start, 1)
        self.assertEqual(interval.end, 10)

    def test_should_initialise_from_chrom_only_region_string(self):
        interval = ChromInterval.from_string("20")
        self.assertEqual(interval.chrom, "20")


class TestChromIntervalListMerging(unittest.TestCase):
    def test_should_merge_identical_intervals(self):
        chrom_interval = interval.ChromInterval("chr1", 1, 10)
        interval_list = [chrom_interval, chrom_interval]

        merged_intervals = merge_intervals(interval_list)
        self.assertListEqual([chrom_interval], merged_intervals)

    def test_should_merge_overlapping_intervals(self):
        interval_list = [
            interval.ChromInterval("chr1", 1, 10),
            interval.ChromInterval("chr1", 8, 15)
        ]

        merged_intervals = merge_intervals(interval_list)
        self.assertListEqual(
            [interval.ChromInterval("chr1", 1, 15)], merged_intervals)

    def test_should_merge_touching_intervals(self):
        interval_list = [
            interval.ChromInterval("chr1", 1, 10),
            interval.ChromInterval("chr1", 10, 15)
        ]

        merged_intervals = merge_intervals(interval_list)
        self.assertListEqual(
            [interval.ChromInterval("chr1", 1, 15)], merged_intervals)

    def test_should_not_merge_non_overlapping_intervals(self):
        interval_list = [
            interval.ChromInterval("chr1", 1, 10),
            interval.ChromInterval("chr1", 11, 15)
        ]

        merged_intervals = merge_intervals(interval_list)
        self.assertListEqual(interval_list, merged_intervals)


class TestReadInterval(unittest.TestCase):
    def test_should_fail_if_string_not_formatted_correctly(self):
        interval_rep = "kasj,a"
        self.assertRaises(Exception, read_interval, interval_rep)

    def test_should_fail_if_start_is_not_castable_to_int(self):
        interval_rep = "xxx-1000"
        self.assertRaises(ValueError, read_interval, interval_rep)

    def test_should_fail_if_end_is_not_castable_to_int(self):
        interval_rep = "0-xxx"
        self.assertRaises(ValueError, read_interval, interval_rep)

    def test_should_fail_if_string_contains_too_many_deliminators(self):
        interval_rep = "0-100-1000"
        self.assertRaises(ValueError, read_interval, interval_rep)

    def test_should_fail_if_end_is_less_or_equal_to_start(self):
        interval_rep = "100-100"
        self.assertRaises(EchidnaException, read_interval, interval_rep)

    def test_should_parse_correctly_formatted_string(self):
        interval_rep = "10-100"
        self.assertEqual(
            read_interval(interval_rep),
            interval.Interval(
                10,
                100))

    def test_should_parse_string_representation_of_interval_class(self):
        test_interval = interval.Interval(19, 20)
        self.assertEqual(
            read_interval(
                str(test_interval)), interval.Interval(
                19, 20))


class TestReadChromInterval(unittest.TestCase):
    def test_should_read_chrom_represenation_correctly(self):
        chrom_representation = "20"

        test_interval = read_chrom_interval(chrom_representation)

        self.assertEqual(
            test_interval,
            interval.ChromInterval(chrom_representation))

    def test_should_fail_if_missing_chrom(self):
        region_string = "0-100"
        self.assertRaises(ValueError, read_chrom_interval, region_string)

    def test_should_read_full_region_string_correctly(self):
        region_string = "20:0-100"
        test_chrom_interval = read_chrom_interval(region_string)
        self.assertEqual(
            test_chrom_interval,
            interval.ChromInterval(
                "20",
                0,
                100))

    def test_should_fail_if_string_contains_too_many_deliminators(self):
        region_string = "20:0-100:20"
        self.assertRaises(ValueError, read_chrom_interval, region_string)


class TestReadChromIntervals(unittest.TestCase):
    def test_should_read_list_of_comma_separated_chrom_intervals_into_list_of_correct_size(self):
        chrom_represenation = ",".join([str(_) for _ in range(0, 2)])
        intervals = list(read_chrom_intervals(chrom_represenation))

        self.assertEqual(len(intervals), 2)


class TestChromIntervalLessThan(unittest.TestCase):
    def test_should_fail_if_chromosome_is_not_standarised(self):
        test_interval_1 = ChromInterval("BAA", 0, 100000)
        test_interval_2 = ChromInterval("20", 0, 100000)
        self.assertTrue(test_interval_2 < test_interval_1)

    def test_should_return_true_if_chrom_is_smaller_in_standard_chrom_ordering(self):
        test_interval_1 = ChromInterval("1", 1000000, 2000000)
        test_interval_2 = ChromInterval("2", 0, 1)
        self.assertTrue(test_interval_1 < test_interval_2)

    def test_should_return_false_if_chrom_is_greater_in_standard_chrom_ordering(self):
        test_interval_1 = ChromInterval("2", 0, 1)
        test_interval_2 = ChromInterval("1", 1000000, 2000000)
        self.assertFalse(test_interval_1 < test_interval_2)

    def test_should_return_true_if_chroms_are_equal_and_interval_compares_less(self):
        test_interval_1 = ChromInterval("20", 0, 1)
        test_interval_2 = ChromInterval("20", 10000000, 20000000)
        self.assertTrue(test_interval_1 < test_interval_2)

    def test_should_return_false_if_chroms_are_equal_and_interval_does_not_compares_less(self):
        test_interval_1 = ChromInterval("20", 10000000, 20000000)
        test_interval_2 = ChromInterval("20", 0, 1)
        self.assertFalse(test_interval_1 < test_interval_2)


class TestChromIntervalGreaterThan(unittest.TestCase):
    def test_should_sort_chromosome_is_not_standarised(self):
        test_interval_1 = ChromInterval("BAA", 0, 100000)
        test_interval_2 = ChromInterval("20", 0, 100000)
        self.assertTrue(test_interval_1 > test_interval_2)

    def test_should_return_false_if_chrom_is_smaller_in_standard_chrom_ordering(self):
        test_interval_1 = ChromInterval("1", 1000000, 2000000)
        test_interval_2 = ChromInterval("2", 0, 1)
        self.assertFalse(test_interval_1 > test_interval_2)

    def test_should_return_true_if_chrom_is_greater_in_standard_chrom_ordering(self):
        test_interval_1 = ChromInterval("2", 0, 1)
        test_interval_2 = ChromInterval("1", 1000000, 2000000)
        self.assertTrue(test_interval_1 > test_interval_2)

    def test_should_return_false_if_chroms_are_equal_and_interval_compares_less(self):
        test_interval_1 = ChromInterval("20", 0, 1)
        test_interval_2 = ChromInterval("20", 10000000, 20000000)
        self.assertFalse(test_interval_1 > test_interval_2)

    def test_should_return_true_if_chroms_are_equal_and_interval_does_not_compares_less(self):
        test_interval_1 = ChromInterval("20", 10000000, 20000000)
        test_interval_2 = ChromInterval("20", 0, 1)
        self.assertTrue(test_interval_1 > test_interval_2)
