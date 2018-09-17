# All content Copyright (C) 2018 Genomics plc

from wecall.utils.interval import Interval


class IntervalRegionIterator(object):

    def __init__(self, iterator, make_interval):
        self.__make_interval = make_interval
        self.__iterator = iterator
        self.__item = None
        self.__interval = None
        self.__previous_results = []
        self.__previous_interval = None

    def __next(self):
        self.__item = next(self.__iterator)
        self.__interval = self.__make_interval(self.__item)

    def __call__(self, interval):
        # validation
        if self.__previous_interval is not None and interval.start < self.__previous_interval.start:
            raise Exception()
        self.__previous_interval = interval

        # pre-process cached reads
        next_results = []
        for new_interval, item in self.__previous_results:
            if new_interval.fast_overlap(interval):
                next_results.append((new_interval, item))

        # add reads to cache
        try:
            if self.__item is None:
                self.__next()
            while True:
                if self.__interval.fast_overlap(interval):
                    next_results.append((self.__interval, self.__item))
                elif self.__interval.end <= interval.start:
                    pass
                else:
                    break
                previous_start = self.__interval.start
                self.__next()
                assert self.__interval.start >= previous_start
        except StopIteration:
            pass

        self.__previous_results = next_results
        return [item for new_interval, item in self.__previous_results]


class BAMRegionIterator(IntervalRegionIterator):

    def __init__(self, fetch_iterator):
        IntervalRegionIterator.__init__(
            self,
            fetch_iterator,
            lambda read: Interval(read.pos, read.aend)
        )
