# All content Copyright (C) 2018 Genomics plc
import ast
import re


class TimingDataItem(object):

    def __init__(self, timing_type, length, length_units, metadata):
        self.timing_type = timing_type
        self.length = length
        self.length_units = length_units
        self.metadata = metadata

    def __eq__(self, other):
        return all((
            self.timing_type == other.timing_type,
            self.length == other.length,
            self.length_units == other.length_units,
            self.metadata == other.metadata,
        ))

    def __ne__(self, other):
        return not self.__eq__(other)


def log_timing_parser(fp):

    log_timing_regex = re.compile(
        "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} -- TIMING -- (?P<payload>.*)\s*$"
    )
    timing_log_lines = []
    for line in fp:
        match = log_timing_regex.match(line)
        if match:
            timing_log_lines.append(match.group('payload'))

    timing_message_regex = re.compile(
        "^(?P<type>\w+)\s+(?P<length>\d+)(?P<units>\w+):\s*"
        "(?P<metadata>(?:[a-zA-z][-_0-9a-zA-Z]*=\"(?:(?:\\\")|[^\"])*\";\s*)*)$")
    errors = []
    timing_data = []
    for timing_message in timing_log_lines:
        match = timing_message_regex.match(timing_message)
        if match:
            metadata_item_regex = re.compile(
                '^(?P<key>[a-zA-z][-_0-9a-zA-Z]*)=(?P<value>\"(?:(?:\\\")|[^\"])*\");$')
            metadata_item_match = metadata_item_regex.match(
                match.group('metadata'))

            timing_data.append(TimingDataItem(
                match.group('type'),
                int(match.group('length')),
                match.group('units'),
                {(metadata_item_match.group('key')): ast.literal_eval(metadata_item_match.group('value'))}
            ))
        else:
            errors.append('failed to parse {!r}'.format(timing_message))
    if errors:
        raise Exception(
            'failed to parse log timings:\n{!r}'.format(
                '\n'.join(errors)))

    return timing_data
