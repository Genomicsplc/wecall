# All content Copyright (C) 2018 Genomics plc
# -*- coding: utf8 -*-

import logging
from math import factorial

from wecall.common.exceptions import weCallException
from wecall.vcfutils.stringutils import from_vcf_str

logger = logging.getLogger(__name__)

# Module level constants
UNKNOWN = '.'


def n_choose_k(n, k):
    return int(factorial(n) / factorial(k) / factorial(n - k))


def _parse_flag(value):
    """
    Parses a 'flag' info field.  If flag is used as a
    proper flag the value is None and it is assumed that
    that means True.  Missing flag is unclear and hence not parsed.
    """
    if value == UNKNOWN:
        return None
    else:
        if value is None:
            return True
        if isinstance(value, bool):
            return value
        value = value.upper()
        if value in {'1', 'YES', 'TRUE'}:
            return True
        elif value in {'0', 'NO', 'FALSE'}:
            return False
        else:
            # For strict VCF parsing configure parser to throw on log warnings.
            # TODO: Work out how to configure logger to do this.
            logging.warning("Invalid flag {}".format(value))
            raise weCallException("Invalid flag {}".format(value))


def make_optional(cls, default):
    def optional(value):
        if value is '.':
            return None
        else:
            try:
                return cls(value)
            except ValueError as ex:
                logger.warn(
                    '{} - using {} as default value'.format(ex, default))
                return default
    return optional


DATATYPE_MAPPING = {
    "Integer": make_optional(int, 0),
    "Float": make_optional(float, 0.0),
    "Flag": _parse_flag,
    "Character": make_optional(str, ''),
    "String": make_optional(str, ''),
}


def parse_cardinality_A(values, number_alts, number_field, number_genotypes):
    if values and len(values) != number_alts:
        logger.warn(
            "Incorrect num ALT values {} instead of {}".format(
                len(values), number_alts))

    # If data missing for Allele in VCF mark it as unknown.
    extracted_data = [[value] for value in values]
    extracted_data.extend([UNKNOWN] * (number_alts - len(extracted_data)))
    return extracted_data


def parse_cardinality_R(values, number_alts, number_field, number_genotypes):
    assert len(values) == number_alts + 1
    return [[values[0], value] for value in values[1:]]


def parse_cardinality_G(values, number_alts, number_field, number_genotypes):
    assert len(values) == number_genotypes
    return [[value] for value in values]


def parse_unknown_cardinality(
        values,
        number_alts,
        number_field,
        number_genotypes):
    return [values for _ in range(number_alts)]


def parse_cardinality(values, number_alts, number_field, number_genotypes):
    if number_field == '0':
        assert values == [True]
    else:
        number = from_vcf_str(number_field, int)
        if len(values) != number:
            logger.debug(
                "Incorrect num values {} instead of {}".format(
                    len(values), number))

    return parse_unknown_cardinality(
        values, number_alts, number_field, number_genotypes)


CARDINALITY_MAPPING = {
    'A': parse_cardinality_A,
    'R': parse_cardinality_R,
    'G': parse_cardinality_G,
    UNKNOWN: parse_unknown_cardinality
}


def pad_vcf_data(parsed_data, expected_number_items):
    if len(parsed_data) != expected_number_items:
        logger.warn(
            'expected {} items in {!r}'.format(
                expected_number_items,
                parsed_data))
    parsed_data.extend(
        [None] for index in range(
            expected_number_items -
            len(parsed_data)))
    return parsed_data


def make_split_info_alt_func(cardinality):
    if cardinality == 'A':
        def split_alts(data, n_alts):
            return pad_vcf_data([[value] for value in data], n_alts)
    elif cardinality == 'R':
        def split_alts(data, n_alts):
            ref = data[0]
            return pad_vcf_data([[ref, alt] for alt in data[1:]], n_alts)
    else:
        def split_alts(data, n_alts):
            return [[item for item in data] for _ in range(n_alts)]
    return split_alts


#  This is to work out the position of genotype likelihood, see VCF spec
def choose_gl_position_for_diploid(first_allele_pos, second_allele_pos):
    return int(second_allele_pos * (second_allele_pos + 1) / 2) + \
        first_allele_pos


def make_split_sample_alt_func(cardinality, parser):
    if cardinality == 'A':
        def split_alts(data, n_alts, gt):
            return pad_vcf_data([[parser(value)] for value in data], n_alts)
    elif cardinality == 'R':
        def split_alts(data, n_alts, gt):
            ref = parser(data[0])
            return pad_vcf_data([[ref, parser(alt)]
                                 for alt in data[1:]], n_alts)
    elif cardinality == 'G':
        def split_alts(data, n_alts, gt):
            if gt is None:
                logger.warn('Unknown ploidy when parsing genotype likelihood')
                return [data for _ in range(n_alts)]

            if len(gt) not in {1, 2}:
                logger.warn(
                    "Unable to handle ploidy other than haploid or diploid.")
                return [data for _ in range(n_alts)]

            expected_number_of_genotype_likelihoods = n_choose_k(
                len(gt) + n_alts, n_alts)

            if len(data) == expected_number_of_genotype_likelihoods:
                if len(gt) == 1:
                    return [[parser(data[0]), parser(data[1 + index])]
                            for index in range(n_alts)]
                elif len(gt) == 2:
                    return [
                        [
                            parser(data[0]),
                            parser(data[choose_gl_position_for_diploid(0, 1 + index)]),
                            parser(data[choose_gl_position_for_diploid(1 + index, 1 + index)])
                        ]
                        for index in range(int(n_alts))
                    ]
            else:
                logger.warn(
                    "Incorrect number of values 'G' cardinality, expected {}, got {}".format(
                        expected_number_of_genotype_likelihoods, len(data)))
                if len(gt) == 1:
                    return [[None, None] for _ in range(n_alts)]
                elif len(gt) == 2:
                    return [[None, None, None] for _ in range(n_alts)]
    else:
        def split_alts(data, n_alts, gt):
            return [[parser(item) for item in data] for _ in range(n_alts)]
    return split_alts


class InfoMetadata(object):

    def __init__(
            self,
            number,
            data_type,
            description,
            source=None,
            version=None):
        self.number = number
        self.data_type = data_type
        self.description = description
        self.source = source
        self.version = version
        self.parser = FieldMetadata(
            None, number, data_type, description, source, version)
        self.split_alts = make_split_info_alt_func(number)

    def __repr__(self):
        data_items = (
            "number={!r}".format(self.number),
            "data_type={!r}".format(self.data_type),
            "description={!r}".format(self.description),
            "source={!r}".format(self.source),
            "version={!r}".format(self.version),
        )
        return "<{!s}: {!s}>".format(
            type(self).__name__, ", ".join(data_items))

    def __eq__(self, other):
        return all((
            self.number == other.number,
            self.data_type == other.data_type,
            self.description == other.description,
            self.source == other.source,
            self.version == other.version,
        ))


class SampleMetadata(object):

    def __init__(self, number, data_type, description):
        self.number = number
        self.data_type = data_type
        self.description = description
        self.parser = FieldMetadata(None, number, data_type, description)
        self.split_alts = make_split_sample_alt_func(
            number, self.parser.parse_func)

    def __repr__(self):
        data_items = (
            "number={!r}".format(self.number),
            "data_type={!r}".format(self.data_type),
            "description={!r}".format(self.description),
        )
        return "<{!s}: {!s}>".format(
            type(self).__name__, ", ".join(data_items))

    def __eq__(self, other):
        return all((
            self.number == other.number,
            self.data_type == other.data_type,
            self.description == other.description,
        ))


class FilterMetadata(object):

    def __init__(self, description):
        self.description = description

    def __eq__(self, other):
        return all((
            self.description == other.description,
        ))

    def __repr__(self):
        return '<{cls}: {val}>'.format(
            cls=type(self).__name__,
            val='{key}={value!r}'.format(
                key='description',
                value=self.description))


class ContigMetadata(object):

    def __init__(self, length=None):
        self.length = length

    def __eq__(self, other):
        return all((
            self.length == other.length,
        ))

    def __repr__(self):
        return '<{cls}: {val}>'.format(
            cls=type(self).__name__,
            val='{key}={value!r}'.format(
                key='length',
                value=self.length))


class AdapterMetadata(object):

    def __init__(self, adapter, hash, date):
        self.adapter = adapter
        self.hash = hash
        self.date = date

    def __eq__(self, other):
        return all((
            self.adapter == other.adapter,
            self.hash == other.hash,
            self.date == other.date,
        ))

    def __repr__(self):
        return '<{cls}: {val}>'.format(
            cls=type(self).__name__, val='; '.join((
                '{key}={value!r}'.format(key=key, value=value) for key, value in [
                    ('adapter', self.adapter),
                    ('hash', self.hash),
                    ('date', self.date),
                ]
            ))
        )


class FieldMetadata(object):
    """
    A class corresponding to either "##FORMAT" or "##INFO"
    line in the vcf header
    """

    # TODO: split into type-specific parsers, improving speed by simplifying
    # code paths

    def __init__(
            self,
            name,
            number_field,
            data_type_field,
            description,
            source=None,
            version=None):
        self.name = name
        self.number_field = number_field
        self.cardinality_func = CARDINALITY_MAPPING.get(
            number_field, parse_cardinality)
        self.data_type_field = data_type_field
        self.parse_func = DATATYPE_MAPPING[data_type_field]
        self.description = description
        self.source = source
        self.version = version

    def __eq__(self, other):
        return (
            self.name == other.name and
            self.number_field == other.number_field and
            self.cardinality_func == other.cardinality_func and
            self.data_type_field == other.data_type_field and
            self.parse_func == other.parse_func and
            self.description == other.description and
            self.source == other.source and
            self.version == other.version
        )

    def __str__(self):
        remainder = ""
        if self.source:
            remainder += ",Source=\"{}\"".format(self.source)
        if self.version:
            remainder += ",Version=\"{}\"".format(self.version)
        return "<ID={},Number={},Type={},Description=\"{}\"{}>".format(
            self.name,
            self.number_field,
            self.data_type_field,
            self.description,
            remainder
        )

    def __call__(self, value):
        return self.parse_func(value)

    def extract_data(self, field, number_alts, number_genotypes):
        if field is None:
            return [self.parse_func(field)]
        elif field == '.':
            return []
        else:
            values = []
            for value in field.split(","):
                if value:
                    try:
                        values.append(
                            self.parse_func(value) if value != UNKNOWN else None)
                    # Tolerate basic faults from badly-formatted header/VCF
                    except ValueError as e:
                        logging.warn(e.message)
                        values.append(None)
            return self.cardinality_func(
                values, number_alts, self.number_field, number_genotypes)
