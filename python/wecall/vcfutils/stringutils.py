# All content Copyright (C) 2018 Genomics plc
from wecall.common.exceptions import EchidnaException


def to_vcf_str(primitive_type):
    if primitive_type is None:
        return "."
    elif isinstance(primitive_type, list) or isinstance(primitive_type, tuple):
        return ','.join(map(to_vcf_str, primitive_type))
    elif isinstance(primitive_type, float):
        return "{:g}".format(primitive_type)
    else:
        return str(primitive_type)


def from_vcf_str(vcf_str, desired_type):
    try:
        return desired_type(vcf_str) if vcf_str != "." else None
    except ValueError:
        raise EchidnaException(
            "Cannot cast {} to {!r}".format(
                vcf_str, desired_type))
