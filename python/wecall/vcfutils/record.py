# All content Copyright (C) 2018 Genomics plc
# Represents a single VCF record with only a single ALT.

import copy
import logging
import re
import sys
from collections import OrderedDict
from itertools import repeat

from wecall.common.exceptions import EchidnaException
from wecall.genomics import variant
from wecall.genomics.variant import Variant
from wecall.vcfutils.fieldmetadata import UNKNOWN
from wecall.vcfutils.genotype_call import GenotypeCall
from wecall.vcfutils.info_data import DeferredInfoData, InfoData, DeferredInfoValue
from wecall.vcfutils.sample_data import SampleData, GENOTYPE_KEY, GENOTYPE_PHRED_LIKELIHOODS_KEY, \
    GENOTYPE_LIKELIHOODS_KEY
from wecall.vcfutils.stringutils import to_vcf_str, from_vcf_str


logger = logging.getLogger(__name__)

CHROM_COL = 0
POS_COL = 1
ID_COL = 2
REF_COL = 3
ALT_COL = 4
QUALITY_COL = 5
FILTER_COL = 6
INFO_COL = 7
FORMAT_COL = 8
SAMPLE_COL = 9


def read_records(schema, line):
    """
    Extracts a sequence of `Record` objects from a single line in a VCF file.
    """
    try:
        cols = [l for l in line.strip().split("\t")]
        for item in generate_records(schema, cols):
            yield item
    except EchidnaException:
        raise
    except Exception:
        _, exc, tb = sys.exc_info()
        new_exc = EchidnaException(
            "while reading record from line {!r}: {!s}".format(
                line, exc.message))
        raise new_exc.__class__(new_exc).with_traceback(tb)


def generate_records(schema, cols):
    alts = cols[ALT_COL].split(',')
    vars = [Variant(cols[CHROM_COL], int(cols[POS_COL]) -
                    1, cols[REF_COL], alt) for alt in alts]

    info_data_list = []
    if len(alts) == 1:
        # deferred parsing is simple with a single alt
        info_data_list.append(
            DeferredInfoData(
                schema,
                lambda: defer_parse_info_field(
                    schema,
                    cols[INFO_COL])))

    else:
        # extract and split info data into lists of length n_alts
        split_info_data = OrderedDict()
        for key, value in parse_info_field(cols[INFO_COL]):
            try:
                info_metadata = schema.get_info_data(key)
            except KeyError:
                split_info_data[key] = [
                    DeferredInfoValue(
                        schema, key, value) for index in range(
                        len(alts))]
            else:
                split_info_data[key] = info_metadata.split_alts(
                    value if isinstance(value, list) else value.split(','), n_alts=len(alts)
                )

        # construct InfoData objects from prepared info data
        for index in range(len(alts)):
            info_dict = OrderedDict([
                (key, values[index]) for key, values in list(split_info_data.items())
            ])
            info_data_list.append(InfoData(schema, info_dict))

    try:
        sample_format = cols[FORMAT_COL].split(':')
    except IndexError:
        sample_data_list = repeat(None)
    else:
        # extract sample format
        split_sample_data = {sample_name: sample_field.split(
            ':') for sample_name, sample_field in zip(schema.samples, cols[SAMPLE_COL:])}

        sample_data_list = [
            SampleData(
                cols[FORMAT_COL].split(':'),
                schema.samples) for _ in alts]
        for sample_name, sample_items in list(split_sample_data.items()):
            split_sample_items = {}

            # extract data from sample fields
            gt = None
            for key, item in zip(sample_format, sample_items):
                try:
                    if key == GENOTYPE_KEY:
                        gt = GenotypeCall(item)
                        values = [
                            GenotypeCall(gt.deliminator().join(
                                # Note: default value should be '.', but
                                # downstream tools aren't good enough to use it
                                {None: '.', 0: '0', 1 + index: '1'}.get(gt_index, '0') for gt_index in gt
                            ))
                            for index in range(len(alts))
                        ]
                    elif key == GENOTYPE_LIKELIHOODS_KEY or key == GENOTYPE_PHRED_LIKELIHOODS_KEY:
                        values = schema.get_sample_data(key).split_alts(
                            item.split(','), len(alts), gt)
                    else:
                        values = schema.get_sample_data(key).split_alts(
                            item.split(','), len(alts), None)
                    split_sample_items[key] = values
                except Exception as e:
                    raise type(e)(
                        "Error parsing field {} for sample {}: {}".format(
                            key, sample_name, e))

            # distribute data to each split sample meta-data container
            for index in range(len(alts)):
                sample_data = sample_data_list[index]
                for key, value in list(split_sample_items.items()):
                    sample_data.add_sample_data(sample_name, key, value[index])

    # generate & return record objects
    for var, info_data, sample_data in zip(
            vars, info_data_list, sample_data_list):
        qual = variant_quality_from_vcf(cols[QUALITY_COL])
        ids = variant_ids_from_vcf(cols[ID_COL])
        filts = filters_from_vcf(cols[FILTER_COL])
        yield Record(schema, var, ids, qual, filts, info_data, sample_data, len(alts) > 1)


info_item_regex = re.compile(r'^(?P<key>[a-zA-Z0-9,_+-]+)(?:=(?P<value>.*))?$')


def parse_info_field(field):
    for item in field.split(';'):
        match = info_item_regex.match(item)
        if match:
            key, value = match.group('key'), match.group('value')
            if value is None:
                yield key, [True]
            else:
                yield key, value


def defer_parse_info_field(schema, field):
    for item in field.split(';'):
        match = info_item_regex.match(item)
        if match:
            key, value = match.group('key'), match.group('value')
            if value is None:
                yield key, [True]
            else:
                yield key, DeferredInfoValue(schema, key, value)


def variant_from_vcf(chrom_column, pos_column, ref_column, alt_column):
    return variant.Variant(
        chrom_column,
        from_vcf_str(
            pos_column,
            int) - 1,
        ref_column,
        alt_column)


def variant_quality_from_vcf(quality_string):
    return from_vcf_str(quality_string, float)


def variant_ids_from_vcf(id_column):
    return set() if id_column == UNKNOWN else set(id_column.split(","))


def vcf_id_entry_from_variant_ids(variant_ids):
    return ",".join(sorted(variant_ids)) if variant_ids else UNKNOWN


def filters_from_vcf(filter_column):
    return set() if filter_column == "PASS" else set(filter_column.split(";"))


def vcf_row_from_variant(
        variant,
        variant_ids=set(),
        quality=None,
        filters=set(),
        info_data=None,
        sample_data=None,
):
    columns = [
        variant.chrom,
        to_vcf_str(variant.one_indexed_pos_from),
        vcf_id_entry_from_variant_ids(variant_ids),
        variant.ref,
        variant.alt,
        to_vcf_str(quality),
        ";".join(filters) if filters else "PASS",
        info_data.to_vcf() if info_data else UNKNOWN,
    ]
    if sample_data is not None:
        columns = columns + sample_data.to_vcf_columns()
    return "\t".join(columns)


def vcf_row_from_record(
        record,
        variant_ids=set(),
        quality=None,
        filters=None
):
    actual_filters = filters if filters else record.filters
    filter_field = ";".join(actual_filters) if actual_filters else "PASS"
    fields = list([
        record.variant.chrom,
        to_vcf_str(record.one_indexed_pos_from),
        vcf_id_entry_from_variant_ids(variant_ids if variant_ids else record.ids),
        record.variant.ref,
        record.variant.alt,
        to_vcf_str(quality if quality else record.quality),
        filter_field,
        record.info.to_vcf() if record.info else UNKNOWN,
    ])
    if record.sample_info is not None:
        fields.extend(record.sample_info.to_vcf_columns())
    return "\t".join(
        fields
    )


def split_GT(composite_GT, number_alts):
    composite_call = GenotypeCall(composite_GT)
    individual_GTs = [str(composite_call.get_genotype_call_for_alt(i + 1))
                      for i in range(0, number_alts)]
    return individual_GTs


def split_MNP_variant(var, include_ref_calls=False):
    if var.type != variant.TYPE_MNP:
        yield var
    else:
        for index, (ref, alt) in enumerate(zip(var.ref, var.alt)):
            if ref != alt or include_ref_calls:
                yield variant.Variant(var.chrom, index + var.pos_from, ref, alt)


def split_MNP_record(record):
    if record.type != variant.TYPE_MNP:
        yield record
    else:
        for var in split_MNP_variant(record.variant):
            new_record = copy.deepcopy(record)
            new_record.variant = var
            yield new_record


def common_prefix_length(lhs, rhs):
    offset = 0
    for lhs, rhs in zip(lhs, rhs):
        if lhs == rhs:
            offset += 1
        else:
            break
    return offset


def trimmed_vcf_ref_alt(ref, alt):
    if len(ref) == 0 or len(alt) == 0:
        raise EchidnaException("VCF format requires non-empty ref and alt")
    if ref == alt and len(ref) > 1:
        raise EchidnaException("VCF requires refcalls of length 1")
    if alt == UNKNOWN or ref == UNKNOWN:
        # VCF allows this to indicate unknown data.
        raise EchidnaException("not dealing with monomorphic variants")
    offset, new_ref, new_alt = trimmed_ref_alt(ref, alt)
    start_context, end_context = 0, 0
    if len(ref) != len(alt) or (not new_ref and not new_alt):
        if offset == 0:
            end_context = 1
        else:
            start_context = 1
    result_ref =\
        ref[offset - start_context:offset] +\
        new_ref +\
        ref[offset + len(new_ref):offset + len(new_ref) + end_context]
    result_alt =\
        alt[offset - start_context:offset] +\
        new_alt +\
        alt[offset + len(new_alt):offset + len(new_alt) + end_context]
    return offset - start_context, result_ref, result_alt


def trimmed_ref_alt(ref, alt):
    end_offset = common_prefix_length(reversed(ref), reversed(alt))
    start_offset = common_prefix_length(
        ref[:len(ref) - end_offset], alt[:len(alt) - end_offset])
    new_ref_len = len(ref) - start_offset - end_offset
    new_alt_len = len(alt) - start_offset - end_offset
    assert(new_ref_len >= 0)
    assert(new_alt_len >= 0)
    return start_offset, ref[start_offset:len(
        ref) - end_offset], alt[start_offset:len(alt) - end_offset]


def trim_variant(var):
    start_offset, ref, alt = trimmed_vcf_ref_alt(var.ref, var.alt)
    return variant.Variant(var.chrom, start_offset + var.pos_from, ref, alt)


def trim_record(record):
    record.variant = trim_variant(record.variant)
    return record


class Record(object):
    """
    Class representing a single variant with all its attributes
    """

    __slots__ = (
        'schema',
        'variant',
        'ids',
        'quality',
        'filters',
        'info',
        'sample_info',
        'from_multi_alt')

    def __init__(
            self,
            schema,
            variant,
            variant_id,
            quality,
            filters,
            info,
            sample_info,
            from_multi_alt
    ):
        self.schema = schema
        self.variant = variant
        self.ids = variant_id
        self.quality = quality
        self.filters = filters
        self.info = info
        self.sample_info = sample_info
        self.from_multi_alt = from_multi_alt

    def __hash__(self):
        return hash((self.variant.__repr__(),))

    def __eq__(self, other):
        return (
            self.variant == other.variant and
            self.ids == other.ids and
            self.quality == other.quality and
            self.filters == other.filters and
            self.info == other.info and
            self.sample_info == other.sample_info and
            self.from_multi_alt == other.from_multi_alt
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        data_items = (
            "variant={!r}".format(self.variant),
            "id={!r}".format(self.ids),
            "quality={!r}".format(self.quality),
            "filters={!r}".format(self.filters),
            "info={!r}".format(self.info),
            "sample_info={!r}".format(self.sample_info),
            "from_multi_alt={!r}".format(self.from_multi_alt),
        )
        return "<Record: {!s}>".format(
            ", ".join(data_items)
        )

    @property
    def insert_size(self):
        return self.variant.insert_size

    @property
    def passes_filter(self):
        return len(self.filters) == 0

    @property
    def chrom(self):
        return self.variant.chrom

    @property
    def length(self):
        return self.variant.length

    @property
    def pos_from(self):
        return self.variant.pos_from

    @property
    def pos_to(self):
        return self.variant.pos_to

    @property
    def one_indexed_pos_from(self):
        return self.variant.one_indexed_pos_from

    @property
    def one_indexed_pos_to(self):
        return self.variant.one_indexed_pos_to

    @property
    def ref(self):
        return self.variant.ref

    @property
    def alt(self):
        return self.variant.alt

    @property
    def type(self):
        return self.variant.type

    @property
    def samples(self):
        """
        :return: List with sample names
        """
        return self.schema.samples[:]

    @property
    def genotypes(self):
        """
        :return Dictionary from sample id to a genotype
        """
        return self.sample_info.genotypes()

    def get_one_based_key(self):
        return self.chrom, self.one_indexed_pos_from, self.ref, self.alt
