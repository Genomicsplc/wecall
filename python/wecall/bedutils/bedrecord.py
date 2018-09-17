# All content Copyright (C) 2018 Genomics plc
import functools

from wecall.genomics.variant import Variant
from wecall.utils.interval import ChromInterval


@functools.total_ordering
class BEDRecord(object):

    def __init__(
            self,
            chrom,
            start,
            end,
            label=None,
            score=None,
            strand=None,
            build=None):
        self.chrom_interval = ChromInterval(chrom, start, end)
        self.build = build
        self.label = label
        self.score = score
        self.strand = strand

    def __hash__(self):
        return hash((
            self.chrom_interval,
            self.build,
            self.label,
            self.score,
            self.strand
        ))

    def __repr__(self):
        return "<{}:{}>".format(
            type(self).__name__,
            ", ".join((
                "chrom_interval: {}".format(self.chrom_interval),
                "build: {}".format(self.build),
                "label: {}".format(self.label),
                "score: {}".format(self.score),
                "strand: {}".format(self.strand),
            ))
        )

    def __len__(self):
        return len(self.chrom_interval)

    def __lt__(self, other):
        return self.chrom_interval < other.chrom_interval

    def __eq__(self, other):
        return all((
            self.chrom_interval == other.chrom_interval,
            self.build == other.build,
            self.label == other.label,
            self.score == other.score,
            self.strand == other.strand
        ))

    @property
    def interval(self):
        return self.chrom_interval.interval

    @property
    def chrom(self):
        return self.chrom_interval.chrom

    @property
    def start(self):
        return self.chrom_interval.start

    @property
    def end(self):
        return self.chrom_interval.end

    @property
    def variant(self):
        return Variant(self.chrom, self.start, '.', '.')

    def __str__(self):
        return "\t".join(
            str(r) if r is not None else '.' for r in [
                self.chrom,
                self.start,
                self.end,
                self.label,
                self.score,
                self.strand])
