# All content Copyright (C) 2018 Genomics plc
from collections import OrderedDict
from wecall.vcfutils import fieldmetadata


class Schema(object):

    def __init__(self):
        # TODO: file format is a property of the serialisation, not the data
        self.vcf_format = None
        self.file_metadata = OrderedDict()
        self.samples = []
        self.__info_data = OrderedDict()
        self.__sample_data = OrderedDict()
        self.__filters = OrderedDict()
        self.__contigs = OrderedDict()
        self.__adapters = []

    def __eq__(self, other):
        return all((
            self.file_metadata == other.file_metadata,
            self.samples == other.samples,
            self.__info_data == other.__info_data,
            self.__sample_data == other.__sample_data,
            self.__filters == other.__filters,
            self.__contigs == other.__contigs,
        ))

    def __repr__(self):
        data_items = (
            "file_metadata={!r}".format(self.file_metadata),
            "samples={!r}".format(self.samples),
            "info_data={!r}".format(self.__info_data),
            "sample_data={!r}".format(self.__sample_data),
            "filters={!r}".format(self.__filters),
            "contigs={!r}".format(self.__contigs)
        )
        return "<{}: {!s}>".format(type(self).__name__, ", ".join(data_items))

    # TODO: remove this function
    def set_vcf_format(self, vcf_format):
        self.vcf_format = vcf_format

    # info data:

    def set_info_data(
            self,
            key,
            number,
            data_type,
            description,
            source=None,
            version=None):
        self.__info_data[key] = fieldmetadata.InfoMetadata(
            number, data_type, description, source, version)

    def get_info_data(self, key):
        return self.__info_data[key]

    def iter_info_data(self):
        for key, value in list(self.__info_data.items()):
            yield key, value

    def del_info_data(self, key):
        del self.__info_data[key]

    # sample data:

    def set_sample_data(self, key, number, data_type, description):
        self.__sample_data[key] = fieldmetadata.SampleMetadata(
            number, data_type, description)

    def get_sample_data(self, key):
        return self.__sample_data[key]

    def iter_sample_data(self):
        for key, value in list(self.__sample_data.items()):
            yield key, value

    def del_sample_data(self, key):
        del self.__sample_data[key]

    # filter data:

    def set_filter(self, key, description):
        self.__filters[key] = fieldmetadata.FilterMetadata(description)

    def get_filter(self, key):
        return self.__filters[key]

    def iter_filters(self):
        for key, value in list(self.__filters.items()):
            yield key, value

    def del_filter(self, key):
        del self.__filters[key]

    # contig data:

    def set_contig(self, key, length=None):
        self.__contigs[key] = fieldmetadata.ContigMetadata(length)

    def get_contig(self, key):
        return self.__contigs[key]

    def iter_contigs(self):
        for key, value in list(self.__contigs.items()):
            yield key, value

    def del_contig(self, key):
        del self.__contigs[key]

    # adapters:

    @property
    def from_adapted_vcf(self):
        return len(self.__adapters) > 0

    def set_adapter(self, adapter, hash, date):
        self.__adapters.append(
            fieldmetadata.AdapterMetadata(
                adapter, hash, date))

    def iter_adapters(self):
        for adapter in self.__adapters:
            yield adapter
