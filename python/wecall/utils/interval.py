# All content Copyright (C) 2018 Genomics plc
import functools
from wecall.common.exceptions import weCallException
from wecall.genomics.chromosome import standardise_chromosome, get_chromosome_index, \
    chromosome_comp


@functools.total_ordering
class Interval(object):
    """
    Stores an interval on a contiguous piece of DNA.
    There is a special null type introduced for algebraic completeness.
    """

    # Magic methods

    def __repr__(self):
        return "<{}>".format(
            ", ".join((repr(x) for x in (
                self.start,
                self.end,
            )))
        )

    def __hash__(self):
        return hash((self.start, self.end))

    def __str__(self):
        return "{!s}-{!s}".format(self.start, self.end)

    def __eq__(self, other):
        return all((self.start == other.start, self.end == other.end))

    def __lt__(self, other):
        if self.start == other.start:
            return self.end < other.end
        else:
            return self.start < other.start

    def __init__(self, start=None, end=None):
        """
        By default a null interval is created.
        """
        if start is None or end is None or start > end:
            self.start = self.end = None
        else:
            self.start = start
            self.end = end
        assert self.__is_valid()

    # Private methods

    def __is_valid(self):
        return (
            all((self.start is None, self.end is None)) or
            not any((self.start is None, self.end is None))
        )

    def __compatible(lhs, rhs):
        assert lhs.__is_valid() and rhs.__is_valid()
        return True

    # Properties

    @property
    def is_null(self):
        """
        A `null` set in this context is an empty set which has an indeterminate
        position.
        """
        assert self.__is_valid()
        return all((self.start is None, self.end is None))

    @property
    def is_empty(self):
        """
        An `empty` set in this context is a set of zero length with a well
        defined position.
        """
        assert self.__is_valid()
        return any((self.start is None, self.end is None)) or self.length == 0

    @property
    def length(self):
        assert self.__is_valid()
        if any((self.start is None, self.end is None)):
            return 0
        else:
            return self.end - self.start

    # Public methods

    def toDict(self):
        return {
            "start": self.start,
            "end": self.end,
        }

    @staticmethod
    def fromDict(dct):
        return Interval(dct["start"], dct["end"])

    def overlap_size(lhs, rhs):
        """
        For compatible intervals returns a quantity of overlap, otherwise raise an exception.
        Returns None if intervals don't overlap.
        """
        if not lhs.__compatible(rhs):
            return None
        if lhs.is_null or rhs.is_null:
            return None
        else:
            size = min(lhs.end, rhs.end) - max(lhs.start, rhs.start)
            if size >= 0:
                return size
            else:
                return None

    def __contains__(self, other):
        return all((
            other.start >= self.start,
            other.end <= self.end,
        ))

    def fast_overlap(lhs, rhs):
        return (min(lhs.end, rhs.end) - max(lhs.start, rhs.start)) > 0

    def overlap(lhs, rhs):
        """
        For compatible intervals returns a boolean value, otherwise raise an exception.
        Returns true if intervals overlap.
        """
        size = Interval.overlap_size(lhs, rhs)
        if size is None:
            return False
        else:
            return size > 0

    def overlapOrTouch(lhs, rhs):
        """
        For compatible intervals returns a boolean value, otherwise raise an exception.
        Returns true if intervals overlap or touch.
        """
        size = Interval.overlap_size(lhs, rhs)
        if size is None:
            return False
        else:
            return size >= 0

    def intersection(lhs, rhs):
        """
        Always returns single interval which may be null when it would
        otherwise be not well-defined.
        Not expected to throw for compatible intervals.
        """
        assert lhs.__compatible(rhs)
        if lhs.is_null or rhs.is_null:
            return Interval()
        else:
            return Interval(
                max(lhs.start, rhs.start),
                min(lhs.end, rhs.end)
            )

    def combination(lhs, rhs):
        """
        Returns lhs with as much of rhs as possible.
        """
        assert lhs.__compatible(rhs)
        if not lhs.overlapOrTouch(
                rhs) and not lhs.end == rhs.start and not rhs.end == lhs.start:
            return lhs
        else:
            return Interval(
                min(lhs.start, rhs.start) if lhs.start is not None and rhs.start is not None else None,
                max(lhs.end, rhs.end) if lhs.end is not None and rhs.end is not None else None
            )

    def leftDifference(lhs, rhs):
        """
        Always returns single interval which may be null when it would
        otherwise be not well-defined.
        """
        assert lhs.__compatible(rhs)
        if lhs.is_null or rhs.is_null:
            return Interval()
        else:
            return Interval(
                lhs.start,
                max(min(rhs.start, lhs.end), lhs.start)
            )

    def rightDifference(lhs, rhs):
        """
        Always returns single interval which may be null when it would
        otherwise be not well-defined.
        """
        assert lhs.__compatible(rhs)
        if lhs.is_null or rhs.is_null:
            return Interval()
        else:
            return Interval(
                min(max(rhs.end, lhs.start), lhs.end),
                lhs.end
            )


def interval_coverage(sub_intervals, main_interval):
    """
    Find length within main_interval that is covered by (overlapping) sub-intervals.
    :param sub_intervals:
    :param main_interval:
    :return:
    """
    length = 0
    interval_to_add = None
    sorted_intervals = sorted(sub_intervals)
    for interval in sorted_intervals:
        current_interval = interval
        if current_interval.overlap(main_interval):
            # Trim the interval
            current_interval = current_interval.intersection(main_interval)

            if interval_to_add is None:
                interval_to_add = current_interval

            if interval_to_add.overlap(current_interval):
                interval_to_add = current_interval.combination(interval_to_add)
            else:
                length += interval_to_add.length
                interval_to_add = current_interval

    if interval_to_add is not None:
        length += interval_to_add.length

    return length


@functools.total_ordering
class ChromInterval(object):
    """
    Class to represent a continuous interval of DNA (on one Chromosome!!)
    """

    def __init__(self, chrom, start=None, end=None):
        """
        """
        self.chrom = chrom
        self.interval = Interval(start, end)

    def __hash__(self):
        return hash((self.chrom, self.interval))

    def __repr__(self):
        return "<ChromInterval: chrom={}, interval={}>".format(
            self.chrom, self.interval)

    def __str__(self):
        if self.interval.is_null:
            # Currently can not run Edna on only second half of chromosome.
            return "{!s}".format(self.chrom)
        else:
            return "{!s}:{!s}".format(self.chrom, self.interval)

    def __eq__(self, other):
        return self.chrom == other.chrom and self.interval == other.interval

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.interval < other.interval
        else:
            return chromosome_comp(self.chrom, other.chrom)

    def __gt__(self, other):
        return not self < other and not self == other

    @classmethod
    def from_string(cls, region_string):
        if ':' in region_string:
            chrom, interval = region_string.split(':')
            start, end = [int(x) for x in interval.split('-')]
            return cls(chrom, start, end)
        else:
            chrom = region_string
            return cls(chrom)

    def overlap(lhs, rhs):
        if lhs.chrom != rhs.chrom:
            return False
        elif lhs.interval.is_null:
            return True
        elif rhs.interval.is_null:
            return True
        else:
            return lhs.interval.overlap(rhs.interval)

    def overlapOrTouch(lhs, rhs):
        if lhs.chrom != rhs.chrom:
            return False
        elif lhs.interval.is_null:
            return True
        elif rhs.interval.is_null:
            return True
        else:
            return lhs.interval.overlapOrTouch(rhs.interval)

    @property
    def chrom_index(self):
        return get_chromosome_index(self.chrom)

    @property
    def start(self):
        return self.interval.start

    @property
    def end(self):
        return self.interval.end

    @end.setter
    def end(self, value):
        self.interval.end = value

    def __len__(self):
        return self.length

    @property
    def length(self):
        return self.interval.length

    def intersection(lhs, rhs):
        assert(lhs.chrom == rhs.chrom)
        # interval.is_null = the whole chromosome.
        if lhs.interval.is_null:
            return rhs
        elif rhs.interval.is_null:
            return lhs
        else:
            new_interval = lhs.interval.intersection(rhs.interval)
            return ChromInterval(
                lhs.chrom,
                new_interval.start,
                new_interval.end
            )

    def normalised(self):
        return ChromInterval(
            standardise_chromosome(
                self.chrom),
            self.interval.start,
            self.interval.end)


WHOLE_GENOME = ChromInterval(None)


def read_interval(interval_string):
    start_string, end_string = tuple(interval_string.split("-"))
    start, end = int(start_string), int(end_string)
    if end <= start:
        raise weCallException(
            "Interval {} does not have start < end".format(interval_string))
    return Interval(start, end)


def read_chrom_interval(region_string):
    if all(deliminator not in region_string for deliminator in ':,.-'):
        return ChromInterval(region_string)
    else:
        chrom, interval_string = tuple(region_string.split(":"))
        # constructor of chrom interval needs to be changed.
        tmp_interval = read_interval(interval_string)
        return ChromInterval(chrom, tmp_interval.start, tmp_interval.end)


def read_full_chrom_interval(region_string, reference_genome):
    if all(deliminator not in region_string for deliminator in ':,.-'):
        chrom = region_string
        return ChromInterval(
            chrom, 0, reference_genome.get_chrom_length(chrom))
    else:
        chrom, interval_string = tuple(region_string.split(":"))
        # constructor of chrom interval needs to be changed.
        tmp_interval = read_interval(interval_string)
        return ChromInterval(chrom, tmp_interval.start, tmp_interval.end)


def read_chrom_intervals(region_line):
    for region_string in region_line.split(","):
        yield read_chrom_interval(region_string)


def read_full_chrom_intervals(region_line, reference_genome):
    for region_string in region_line.split(","):
        yield read_full_chrom_interval(region_string, reference_genome)


def merge_intervals(sorted_intervals):
    it = iter(sorted_intervals)
    try:
        merged_intervals = [next(it)]
    except StopIteration:
        return []
    for interval in it:
        if merged_intervals[-1].overlapOrTouch(interval):
            merged_intervals[-1].end = interval.end
        else:
            merged_intervals.append(interval)

    return merged_intervals
