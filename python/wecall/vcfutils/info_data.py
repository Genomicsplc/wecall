# All content Copyright (C) 2018 Genomics plc
import logging
from collections import OrderedDict

from wecall.vcfutils.fieldmetadata import UNKNOWN

logger = logging.getLogger()


class DeferredInfoData(object):

    __slots__ = ('__schema', '__data_generator', '__data')

    def __init__(self, schema, data_generator):
        self.__schema = schema
        self.__data_generator = data_generator
        self.__data = None

    def __repr__(self):
        serialised_data = '(unparsed)' if self.__data is None else ', '.join(
            '{}: {}'.format(key, value) for key, value in list(self.__data.items())
        )
        return '<{}: {}>'.format(type(self).__name__, serialised_data)

    def __require_data_mapping(self):
        if self.__data is None:
            self.__data = OrderedDict(self.__data_generator())

    def __require_all_values(self):
        self.__require_data_mapping()
        for key, value in list(self.__data.items()):
            if isinstance(value, DeferredInfoValue):
                self.__data[key] = value()

    def __len__(self):
        self.__require_data_mapping()
        return len(self.__data)

    def __contains__(self, key):
        self.__require_data_mapping()
        return key in self.__data

    def keys(self):
        self.__require_data_mapping()
        return list(self.__data.keys())

    def values(self):
        self.__require_all_values()
        return list(self.__data.values())

    def items(self):
        self.__require_all_values()
        return list(self.__data.items())

    def __getitem__(self, key):
        self.__require_data_mapping()
        value = self.__data[key]
        if isinstance(value, DeferredInfoValue):
            value = self.__data[key] = value()
        return value

    def __setitem__(self, key, value):
        self.__require_data_mapping()
        self.__data[key] = value

    def __eq__(self, other):
        return dict(list(self.items())) == dict(list(other.items()))

    def to_vcf(self):
        self.__require_all_values()
        if len(self.__data) == 0:
            return UNKNOWN
        else:
            info_strings = []
            for key, value in list(self.__data.items()):
                if value is None:
                    info_strings.append(key)
                else:
                    info_strings.append("{!s}={!s}".format(
                        key, ','.join(map(str, value))))
            return ";".join(info_strings)


class DeferredInfoValue(object):

    __slots__ = ('__schema', '__key', '__value', '__parser')

    def __init__(self, schema, key, value):
        self.__schema = schema
        self.__key = key
        self.__value = value
        try:
            info_data = self.__schema.get_info_data(self.__key)
        except KeyError:
            self.__parser = None
        else:
            self.__parser = info_data.parser

    def __call__(self):
        if self.__parser is None:
            try:
                info_data = self.__schema.get_info_data(self.__key)
            except KeyError:
                logger.warn(
                    'info field {!r} not defined in schema'.format(
                        self.__key))
                self.__schema.set_info_data(
                    self.__key,
                    '.',
                    'String',
                    'Inferred from file content during parsing',
                    'vcfutils',
                    'undefined')
                info_data = self.__schema.get_info_data(self.__key)
            self.__parser = info_data.parser
        try:
            return [self.__parser(item) for item in self.__value.split(',')]
        except Exception as e:
            raise type(e)(
                "Error parsing field {}: {}".format(
                    self.__key, e.message))


class InfoData(object):
    """
    Class that represents all the info fields. Acts as a dictionary INFO_KEY: INFO_VALUE
    """

    __slots__ = ('__schema', '__dict')

    def __init__(self, schema, dict):
        self.__schema = schema
        self.__dict = dict

    def to_vcf(self):
        if len(self.__dict) == 0:
            return UNKNOWN
        else:
            info_strings = []
            for key, value in sorted(self.__dict.items()):
                if value is None:
                    info_strings.append(key)
                else:
                    info_strings.append("{!s}={!s}".format(
                        key, ','.join(map(str, value))))
            return ";".join(info_strings)

    def __setitem__(self, key, value):
        try:
            info_data = self.__schema.get_info_data(key)
        except KeyError:
            raise KeyError(
                "Attempt to read INFO field {!r} which is not defined in the VCF header".format(key))
        self.__dict[key] = value

    def __getitem__(self, key):
        value = self.__dict[key]
        if isinstance(value, DeferredInfoValue):
            value = self.__dict[key] = value()
        return value

    def __contains__(self, key):
        return key in self.__dict

    def __len__(self):
        return len(self.__dict)

    def __eq__(self, other):
        return dict(list(self.items())) == dict(list(other.items()))

    def keys(self):
        return list(self.__dict.keys())

    def values(self):
        return list(self.__dict.values())

    def items(self):
        return list(self.__dict.items())

    def __repr__(self):
        serialised_data = ', '.join(
            ("{!r}: {!r}".format(
                key, value) for key, value in list(
                self.__dict.items())))
        return "<InfoData: {" + serialised_data + "}>"
