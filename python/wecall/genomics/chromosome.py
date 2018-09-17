# All content Copyright (C) 2018 Genomics plc
from wecall.common.exceptions import EchidnaException

# Define a standard ordering for chromosomes. N.B. We always use the naming convention "1",
# "2", "3" etc for internal use, not "chr1", "chr2", "chr3"

CHROMOSOME_LIST = [
    '1',
    '2',
    '3',
    '4',
    '5',
    '6',
    '7',
    '8',
    '9',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '18',
    '19',
    '20',
    '21',
    '22',
    'X',
    'Y',
    'MT']
CHROMOSOME_ORDER = {
    chrom: index for index,
    chrom in enumerate(CHROMOSOME_LIST)}


def chromosome_comp(lhs, rhs):
    p_left = CHROMOSOME_ORDER.get(lhs, len(CHROMOSOME_ORDER))
    p_right = CHROMOSOME_ORDER.get(rhs, len(CHROMOSOME_ORDER))
    if p_left != p_right:
        return p_left < p_right
    else:
        return lhs < rhs


def standardise_chromosome(chrom):
    stripped_chrom = chrom.upper().replace("CHR", "").lstrip('0')
    if stripped_chrom == "M":
        return "MT"
    else:
        return stripped_chrom


def add_chr(chrom):
    return "chr{}".format(chrom)


def get_chromosome_index(chrom):
    try:
        return CHROMOSOME_ORDER[standardise_chromosome(chrom)]
    except KeyError:
        raise EchidnaException("Invalid chromosome {}".format(chrom))
