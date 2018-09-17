# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.read_sequence import HIGH_QUALITY
from wecall.bamutils.sequence_builder import sequence_builder


class SequenceBank(object):
    """
    A container to hold annotated DNA sequences in relation to a reference sequence.
    """

    def __init__(self, reference):
        self.reference = reference
        self._read_sequences_with_coverage = []

    def __getitem__(self, item):
        return self._read_sequences_with_coverage[item]

    def __len__(self):
        return len(self._read_sequences_with_coverage)

    @property
    def chrom(self):
        return self.reference.chrom

    @property
    def variants(self):
        variants = set()
        for sequence in self._read_sequences_with_coverage:
            variants.update(sequence.read_sequence.variants)
        return variants

    def add_sequence(
            self,
            seq_string,
            quality_string=None,
            n_fwd=None,
            n_rev=None,
            mapping_quality=HIGH_QUALITY,
            insert_size=None,
            read_id=None,
            read_flags=None,
            cigar_string=None,
            read_start=None,
            read_mate_start=None
    ):
        self._read_sequences_with_coverage.extend(
            sequence_builder(
                self.reference,
                seq_string,
                quality_string,
                n_fwd,
                n_rev,
                mapping_quality,
                insert_size,
                read_id,
                read_flags,
                cigar_string,
                read_start,
                read_mate_start
            )
        )
        return self

    def build_reads(self, chrom_id, read_tags):
        for read_seq_with_coverage in self._read_sequences_with_coverage:
            for read in read_seq_with_coverage.build_reads(
                    chrom_id, read_tags):
                yield read


class AsciiVariantGenerator(object):

    def __init__(self, reference):
        self.reference = reference

    def get_variants(self, ascii_haplotypes):
        seq_bank = SequenceBank(self.reference)
        for candidate_ascii_haplotype in ascii_haplotypes:
            seq_bank.add_sequence(candidate_ascii_haplotype)
        return seq_bank.variants
