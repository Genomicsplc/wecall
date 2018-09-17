#!/usr/bin/env python
# All content Copyright (C) 2018 Genomics plc
from __future__ import print_function
import sys
import re
import ast
import argparse
import os


print('WARNING: {!r} is obsolete.'.format(os.path.abspath(__file__)))


def main(args):
    exe, template_path, target_path, substitutions = parse_args(args)
    with open(template_path, "r") as template_fp:
        with open(target_path, "w") as target_fp:
            render_template(template_fp, target_fp, substitutions)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('input_template', help='Template to be rendered')
    parser.add_argument('output_target', help='Rendered template')
    parser.add_argument(
        'substitutions',
        nargs="+",
        help="Substitutions in the form +Field=Value")
    exe = args[0]
    parsed_args = parser.parse_args(args[1:])
    substitutions = parse_substitutions(parsed_args.substitutions)
    return exe, parsed_args.input_template, parsed_args.output_target, substitutions


def parse_substitutions(args):
    '''python scripts/renderTemplate a b +a=" " +version=1.2.3 +hehe[] $(ls /) \;g \; +blah=4

    and the example template to go with the above

    {hehe:" "}
    {blah}
    {version}
    '''

    single_value_regex = re.compile(
        '^\\+(?P<field>[a-zA-Z_][a-zA-Z0-9_]*)=(?P<value>.*)$', re.S)

    multiple_value_regex = re.compile(
        '^\\+(?P<field>[a-zA-Z_][a-zA-Z0-9_]*)\\[\\]$', re.S)
    multiple_value_terminator_regex = re.compile('^;$', re.S)

    substitutions = {}
    it = iter(args)
    for token in it:
        single_value_match = single_value_regex.match(token)
        if single_value_match:
            substitutions[single_value_match.group(
                'field')] = single_value_match.group('value')
            continue
        multiple_value_match = multiple_value_regex.match(token)
        if multiple_value_match:
            array = ArraySubstitution()
            substitutions[multiple_value_match.group('field')] = array
            for array_element_token in it:
                if multiple_value_terminator_regex.match(array_element_token):
                    break
                else:
                    array += array_element_token

            continue
        raise Exception(
            'unexpected extra argument: {token!r}'.format(
                token=token))
    return substitutions


class ArraySubstitution(object):

    def __init__(self):
        self.__data = []

    def __repr__(self):
        return('<{name!s}: {data!r}>'.format(name=type(self).__name__, data=self.__data))

    # replacement = '{field:foo}'; format_spec = "foo"
    def __format__(self, format_spec=None):
        maybeMatch = re.match('^(".*")$', format_spec, re.S)
        if maybeMatch:
            separator = ast.literal_eval(maybeMatch.groups()[0])
        else:
            raise Exception(
                "No match while trying to parse separator out of format_spec")

        return separator.join(self.__data)

    def __iadd__(self, element):
        self.__data.append(element)
        return self

    def __eq__(self, other):
        return self.__data == other.__data


def render_template(template_fp, target_fp, substitutions):
    print(template_fp.read().format(**substitutions), file=target_fp)


if __name__ == "__main__":
    main(sys.argv)
