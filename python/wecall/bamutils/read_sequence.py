# All content Copyright (C) 2018 Genomics plc
from pysam import AlignedRead

DUPLICATE = 1024
HIGH_QUALITY = 80
FORWARD_GOOD_READ = 1 + 2 + 32 + 64
REVERSE_GOOD_READ = 1 + 2 + 16 + 128


class ReadSequence(object):

    def __init__(
            self,
            sequence,
            quality,
            mapping_quality,
            insert_size=None,
            read_id=None,
            read_flags=None,
            read_start=None,
            read_mate_start=None
    ):
        self.chrom_id = None
        self.sequence = sequence
        self.quality = quality
        self.mapping_quality = mapping_quality
        self.insert_size = insert_size if insert_size is not None else int(
            2 * self.sequence.length_minus_deletions())
        self.read_flags = read_flags
        self.read_id = read_id if read_id is not None else "no_id"
        self.read_start = read_start if read_start is not None else self.sequence.pos_from

        self.read_mate_start = read_mate_start
        if self.read_mate_start is None:
            read_mate_end = self.read_start + self.insert_size
            self.read_mate_start = read_mate_end - self.sequence.length_minus_deletions()

    @property
    def variants(self):
        return self.sequence.variants

    def build_read(self, read_tags, is_forward):
        read = AlignedRead()

        read.seq = self.sequence.sequence_minus_deletions()
        read.rname = self.chrom_id
        read.pos = self.read_start
        read.mapq = self.mapping_quality
        read.cigarstring = self.sequence.cigar
        read.rnext = self.chrom_id
        read.pnext = self.read_mate_start
        read.tlen = self.insert_size
        read.qual = self.quality.ascii_quality
        read.tags = read_tags

        read.qname = self.read_id

        if self.read_flags is None:
            read.flag = FORWARD_GOOD_READ if is_forward else REVERSE_GOOD_READ
        else:
            read.flag = self.read_flags

        return read


class ReadSequenceWithCoverage(object):

    def __init__(self, read_sequence, n_fwd, n_rev):
        self.read_sequence = read_sequence
        self.n_fwd = n_fwd
        self.n_rev = n_rev

    def build_reads(self, chrom_id, read_tags):
        self.read_sequence.chrom_id = chrom_id
        for fwd_index in range(self.n_fwd):
            yield self.read_sequence.build_read(read_tags, True)
        for rev_index in range(self.n_rev):
            yield self.read_sequence.build_read(read_tags, False)


class ReadPairWithCoverage(object):

    def __init__(self, fwd, rev, n_fwd, n_rev):
        self.__fwd = fwd
        self.__rev = rev
        if n_fwd is not None and n_rev is not None:
            self.n_fwd = n_fwd
            self.n_rev = n_rev
        else:
            self.n_fwd = 1
            self.n_rev = 0

    def build_reads(self, chrom_id, read_tags):
        self.__fwd.chrom_id = chrom_id
        self.__rev.chrom_id = chrom_id
        for fwd_index in range(self.n_fwd):
            first = self.__fwd.build_read(read_tags, True)
            second = self.__rev.build_read(read_tags, True)
            first.mpos = second.pos
            second.mpos = first.pos
            yield first
            yield second
        for rev_index in range(self.n_rev):
            first = self.__fwd.build_read(read_tags, False)
            second = self.__rev.build_read(read_tags, False)
            first.mpos = second.pos
            second.mpos = first.pos
            yield first
            yield second
