# All content Copyright (C) 2018 Genomics plc
from wecall.bamutils.read_sequence import ReadSequence, ReadSequenceWithCoverage, HIGH_QUALITY, ReadPairWithCoverage
from wecall.bamutils.raw_string_sequence import RawStringSequences
from wecall.bamutils.sequence import Sequence
from wecall.bamutils.sequence_position import SequencePosition
from wecall.bamutils.sequence_quality import SequenceQuality
from wecall.common.exceptions import weCallException
from wecall.genomics.reference_chromosome import ReferenceChromosome


def sequence_builder(
        reference,
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
        read_mate_start=None,
):
    quality_string = " " * \
        len(seq_string) if quality_string is None else quality_string
    if not all(
        i is None for i in [
            n_fwd,
            n_rev]) and any(
            i is None for i in [
                n_fwd,
                n_rev]):
        raise weCallException(
            "Invalid combination of forward and reverse reads: n_fwd = {}, n_rev = {} ".format(
                n_fwd, n_rev))

    if len(seq_string) != reference.length_with_deletions():
        raise weCallException(
            "Sequence has to be of the same length as reference. seq_length {}, ref_length {}".format(
                len(seq_string), reference.length_with_deletions()))

    if len(quality_string) != reference.length_with_deletions():
        raise weCallException(
            "Quality string has to be of the same length as reference.")

    ref_pos = reference.pos_from
    current_raw_seq = RawStringSequences()
    sequences = []
    for ref_char, seq_char, qual_char in zip(
            reference.ref_seq, seq_string, quality_string):
        seq_position = SequencePosition(ref_char, seq_char, qual_char)

        if seq_position.is_gap and current_raw_seq.is_ongoing:
            sequences.append(current_raw_seq)
            current_raw_seq = RawStringSequences()
        elif not seq_position.is_gap:
            current_raw_seq.add_position(seq_position, ref_pos)

        ref_pos = seq_position.update_ref_pos(ref_pos)

    if current_raw_seq.is_ongoing:
        sequences.append(current_raw_seq)

    annotated_seqs = []
    if (
        len(sequences) %
        2 == 0 and all(
            (sequences[index].is_forward_seq() for index in range(
            0,
            len(sequences),
            2))) and all(
                (sequences[index].is_reverse_seq() for index in range(
                    1,
                    len(sequences),
                    2)))):
        # sequence of read pairs
        pairs = list(zip(
            (sequences[index] for index in range(0, len(sequences), 2)),
            (sequences[index] for index in range(1, len(sequences), 2))
        ))
        for fwd, rev in pairs:
            annotated_seqs.extend(
                build_annotated_pair(
                    fwd,
                    rev,
                    n_fwd,
                    n_rev,
                    mapping_quality,
                    insert_size,
                    read_id,
                    read_flags,
                    cigar_string,
                    read_start,
                    read_mate_start))
    else:
        # unpaired reads
        for seq in sequences:
            annotated_seqs.extend(
                seq.build_annotated_seq(
                    n_fwd,
                    n_rev,
                    mapping_quality,
                    insert_size,
                    read_id,
                    read_flags,
                    cigar_string,
                    read_start,
                    read_mate_start))
    return annotated_seqs


def build_annotated_pair(
        fwd,
        rev,
        n_fwd,
        n_rev,
        mapping_quality,
        insert_size,
        read_id,
        read_flags,
        cigar_string,
        read_start,
        read_mate_start):
    fwd_reference = ReferenceChromosome(fwd.reference_string, fwd.pos_from)
    rev_reference = ReferenceChromosome(rev.reference_string, rev.pos_from)
    fwd_sequence = Sequence(
        fwd_reference, fwd.sequence_string.replace(
            ",", ".").upper(), cigar_string)
    rev_sequence = Sequence(
        rev_reference, rev.sequence_string.replace(
            ",", ".").upper(), cigar_string)
    fwd_quality = SequenceQuality(fwd.quality_string)
    rev_quality = SequenceQuality(rev.quality_string)

    fwd_read_sequence = ReadSequence(
        fwd_sequence,
        fwd_quality,
        mapping_quality,
        insert_size,
        read_id,
        read_flags,
        read_start,
        read_mate_start)
    rev_read_sequence = ReadSequence(
        rev_sequence,
        rev_quality,
        mapping_quality,
        insert_size,
        read_id,
        read_flags,
        read_start,
        read_mate_start)
    return [
        ReadPairWithCoverage(
            fwd_read_sequence,
            rev_read_sequence,
            n_fwd,
            n_rev)]
