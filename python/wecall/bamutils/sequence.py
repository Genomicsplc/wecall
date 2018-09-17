# All content Copyright (C) 2018 Genomics plc
import re
from wecall.bamutils.cigar import Cigar
from wecall.bamutils.sequence_position import DELETED_BASE, MATCHING_BASE
from wecall.common.exceptions import EchidnaException
from wecall.genomics.variant import Variant, TYPE_DEL, TYPE_INS, TYPE_SNP, TYPE_TO_STR, TYPE_REF


class Sequence(object):
    """
    Object to represent a sequence of DNA.
    """

    def __init__(self, reference, seq_string, cigar_string=None):
        self.__validate_input(reference, seq_string)

        self._reference = reference
        self._seq = seq_string
        self.pos_from = self._reference.pos_from
        self.pos_to = self._reference.pos_to

        self.variants = self._get_variants()
        self.cigar = cigar_string or str(self._get_cigar())

    @property
    def cigar_string(self):
        return str(self.cigar)

    @property
    def raw_sequence(self):
        return self._seq

    def length_minus_deletions(self):
        return len(self.sequence_minus_deletions())

    def length_with_deletions(self):
        return len(self._seq)

    def sequence_with_deletions(self):
        seq = ""
        for ref_char, alt_char in zip(self._reference.ref_seq, self._seq):
            if alt_char == MATCHING_BASE:
                seq += ref_char
            else:
                seq += alt_char
        return seq

    def sequence_minus_deletions(self):
        return self.sequence_with_deletions().replace(DELETED_BASE, "")

    def _get_variants(self):
        variants = set()

        ref_index = self.pos_from - 1
        current_variant = None

        for ref_char, alt_char in zip(self._reference.ref_seq, self._seq):
            if ref_char != DELETED_BASE:
                ref_index += 1

            if ref_char == DELETED_BASE and alt_char == MATCHING_BASE:
                raise EchidnaException(
                    "Invalid sequence at ref position {}".format(ref_index))
            elif ref_char == DELETED_BASE and alt_char == DELETED_BASE:
                continue
            elif alt_char == MATCHING_BASE:
                current_variant = self.__add_variant_to_set(
                    current_variant, None, variants)
                continue

            if ref_char == DELETED_BASE:
                # insertion
                var_pos = ref_index
                var_ref = self._reference[var_pos]
                var_alt = var_ref + alt_char
            elif alt_char == DELETED_BASE:
                # deletion
                var_pos = ref_index - 1
                var_ref = self._reference[var_pos] + ref_char
                var_alt = self._reference[var_pos]
            else:
                # SNP
                var_pos = ref_index
                var_ref = ref_char
                var_alt = alt_char

            new_variant = Variant(
                self._reference.chrom, var_pos, var_ref, var_alt)
            current_variant = self.__add_variant_to_set(
                current_variant, new_variant, variants)

        self.__add_variant_to_set(current_variant, None, variants)
        variants = self.__remove_deletions_from_edges(variants)
        return variants

    @staticmethod
    def __add_variant_to_set(current_variant, new_variant, var_set):
        var_1, current_variant = Sequence.__potentially_merge_adjacent_variants(
            current_variant, new_variant)
        if var_1 is not None and var_1.type != TYPE_REF:
            var_set.add(var_1)
        return current_variant

    def __remove_deletions_from_edges(self, variants):
        filtered_variants = set()
        for var in variants:
            if var.type == TYPE_DEL and (
                    var.pos_from == -1 or var.pos_to == self.pos_to):
                continue
            else:
                filtered_variants.add(var)

        return filtered_variants

    @staticmethod
    def __potentially_merge_adjacent_variants(var_1, var_2):
        if var_1 is None or var_2 is None or var_1.type != var_2.type:
            return var_1, var_2
        else:
            if var_1.type == TYPE_SNP or var_1.type == TYPE_REF:
                return var_1, var_2
            elif var_1.type == TYPE_DEL:
                merged_variant = Variant(
                    var_1.chrom, var_1.pos_from, var_1.ref + var_2.ref[-1], var_1.alt)
                return None, merged_variant
            elif var_1.type == TYPE_INS:
                merged_variant = Variant(
                    var_1.chrom, var_1.pos_from, var_1.ref, var_1.alt + var_2.alt[-1])
                return None, merged_variant
            else:
                raise EchidnaException(
                    "Unexpected variant type: " + TYPE_TO_STR[var_1.type])

    def _get_cigar(self):
        cigar = Cigar()
        for ref_char, alt_char in zip(self._reference.ref_seq, self._seq):
            if alt_char == MATCHING_BASE:
                cigar += Cigar([(Cigar.MATCH, 1)])
            elif alt_char != DELETED_BASE and ref_char != DELETED_BASE:
                # SNP
                cigar += Cigar([(Cigar.MATCH, 1)])
            elif alt_char == DELETED_BASE and ref_char != DELETED_BASE:
                cigar += Cigar([(Cigar.DELETION, 1)])
            elif alt_char != DELETED_BASE and ref_char == DELETED_BASE:
                cigar += Cigar([(Cigar.INSERTION, 1)])
        return cigar

    @staticmethod
    def __validate_input(ref, seq):
        if len(seq) != ref.length_with_deletions():
            raise EchidnaException(
                "Sequence has to be of the same length as reference.")
        if not re.match(r'^[ACGTURYKMSWBDHVN\*\.]*\Z', seq):
            raise EchidnaException(
                "Illegal character in sequence {!r}".format(seq))
