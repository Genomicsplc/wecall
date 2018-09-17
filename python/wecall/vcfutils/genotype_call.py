# All content Copyright (C) 2018 Genomics plc
import copy
from collections import Counter
from fractions import gcd

from wecall.vcfutils.stringutils import from_vcf_str, to_vcf_str


class GenotypeCall(object):

    def __init__(self, genotype_string, normalize=True):
        self.phased = "|" in genotype_string
        self.__alleles = [
            from_vcf_str(
                x, int) for x in genotype_string.split(
                self.deliminator())]
        if normalize is True:
            self.normalise()

    @property
    def alleles(self):
        return self.__alleles

    @property
    def normalized_allele_count(self):
        counter = Counter()
        counter[0] = 0
        for allele in self.__alleles:
            counter[allele] += 1
        try:
            counter[0] += counter[None]
            del counter[None]
        except KeyError:
            pass
        gcd_value = 0
        for value in list(counter.values()):
            gcd_value = gcd(value, gcd_value)
        if gcd_value > 1:
            for key, value in list(counter.items()):
                counter[key] = value / gcd_value
        return tuple(sorted(counter.values()))

    def deliminator(self):
        return "|" if self.phased else "/"

    def __hash__(self):
        return hash(str(self))

    def __len__(self):
        return len(self.__alleles)

    def __getitem__(self, item):
        return self.__alleles[item]

    def __setitem__(self, key, value):
        self.__alleles[key] = value
        self.normalise()

    def __eq__(self, other):
        if not self.phased or not other.phased:
            return sorted(
                self.__alleles,
                key=lambda x: x if x is not None else -
                1) == sorted(
                other.__alleles,
                key=lambda x: x if x is not None else -
                1)
        else:
            return self.__alleles == other.__alleles

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return self.deliminator().join(to_vcf_str(s) for s in self.__alleles)

    def __repr__(self):
        return "GenotypeCall<Phased={}, Alleles={}>".format(
            self.phased, self.__alleles)

    def is_haploid(self):
        return len(self.__alleles) == 1

    def is_diploid(self):
        return len(self.__alleles) == 2

    def is_unknown(self):
        return all((allele is None for allele in self.__alleles))

    def is_called(self):
        return set(self.__alleles) - {0, None} != set()

    def is_heterozygous(self):
        return not self.is_unknown() and len(set(self.__alleles)) == 2

    def is_homozygous(self):
        return not self.is_unknown() and len(set(self.__alleles)) == 1

    def is_homozygous_alt(self):
        return not self.is_unknown() and self.is_homozygous(
        ) and self.__alleles[0] != 0

    def is_homozygous_ref(self):
        return not self.is_unknown() and self.is_homozygous() and not self.is_homozygous_alt()

    def normalise(self):
        if not self.phased:
            self.__alleles.sort(key=lambda x: x if x is not None else -1)

    def get_genotype_call_for_alt(self, alt_number):
        new_genotype = copy.deepcopy(self)
        for index, allele in enumerate(self.__alleles):
            if allele == alt_number:
                new_genotype[index] = 1
            elif allele is not None:
                new_genotype[index] = 0

        return new_genotype


def merge_genotype_calls(genotype_call1, genotype_call2):
    if len(genotype_call1) != 2 or len(genotype_call2) != 2:
        raise Exception("Cannot merge non-diploid genotype calls.")
    elif genotype_call1.is_homozygous_ref():
        return copy.deepcopy(genotype_call2)
    elif genotype_call2.is_homozygous_ref():
        return copy.deepcopy(genotype_call1)
    elif genotype_call1.is_homozygous_alt() and genotype_call2.is_homozygous_alt():
        return GenotypeCall("1/1")
    elif genotype_call1.is_homozygous_alt() or genotype_call2.is_homozygous_alt():
        raise Exception("Cannot merge any homozygous alt genotype calls.")
    elif (genotype_call1.phased and genotype_call2.phased and genotype_call1.is_heterozygous() and
          genotype_call1.is_heterozygous() and genotype_call1 == genotype_call2):
        return copy.deepcopy(genotype_call1)
    else:
        return GenotypeCall("1/1")
