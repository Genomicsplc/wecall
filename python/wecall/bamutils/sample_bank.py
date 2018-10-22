# All content Copyright (C) 2018 Genomics plc
from collections import OrderedDict
from wecall.bamutils.read_sequence import HIGH_QUALITY
from wecall.bamutils.sequence_bank import SequenceBank
from wecall.bamutils.sequence_quality import SequenceQuality
from wecall.common.exceptions import weCallException
from wecall.genomics.reference_chromosome import DEFAULT_CHROM, ReferenceChromosome


class SampleBank(object):

    def __init__(self, reference_string, pos_from=0, chrom=DEFAULT_CHROM):
        self.reference = ReferenceChromosome(reference_string, pos_from, chrom)
        self.__samples = OrderedDict()

    def __getitem__(self, sample_name):
        return self.__samples[sample_name]

    def __len__(self):
        return len(self.__samples)

    def get(self, sample_name, default=None):
        try:
            return self[sample_name]
        except KeyError:
            return default

    @property
    def sample_names(self):
        return list(self.__samples.keys())

    @property
    def variants(self):
        vars = set()
        for seq_bank in list(self.__samples.values()):
            vars.update(seq_bank.variants)

        return vars

    def add_sample_name(self, sample_name):
        if sample_name in self.__samples:
            raise weCallException(
                "Sample {} already exists in the SampleBank.".format(sample_name))
        sequence_bank = SequenceBank(self.reference)
        self.__samples[sample_name] = sequence_bank
        return sequence_bank

    def add_sample_with_seqs_and_quals(
            self,
            sample_name,
            seq_qual_list,
            n_fwd=None,
            n_rev=None,
            mapping_quality=HIGH_QUALITY
    ):
        self.add_sample_name(sample_name)

        current_seq = None
        for look_ahead in seq_qual_list:
            if current_seq is None:
                current_seq = look_ahead
                continue

            if SequenceQuality.is_valid_qual(look_ahead):
                self.__samples[sample_name].add_sequence(
                    current_seq,
                    quality_string=look_ahead,
                    n_fwd=n_fwd,
                    n_rev=n_rev,
                    mapping_quality=mapping_quality
                )
                current_seq = None
            else:
                self.__samples[sample_name].add_sequence(
                    current_seq, n_fwd=n_fwd, n_rev=n_rev, mapping_quality=mapping_quality)
                current_seq = look_ahead

        if current_seq is not None:
            self.__samples[sample_name].add_sequence(
                current_seq, n_fwd=n_fwd, n_rev=n_rev, mapping_quality=mapping_quality)
