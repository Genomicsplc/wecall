# All content Copyright (C) 2018 Genomics plc
import re
import sys


def latex_encode(unencoded):
    return unencoded \
        .replace('\\', '\\textbackslash{}') \
        .replace('_', '\\_') \
        .replace('#', '\\#') \
        .replace('$', '\\$') \
        .replace('%', '\\%') \
        .replace('&', '\\&') \
        .replace('{', '\\{') \
        .replace('}', '\\}') \
        .replace('\\textbackslash\\{\\}', '\\textbackslash{}') \
        .replace('^', '\\textasciicircum{}') \
        .replace('~', '\\textasciitilde{}') \
        .replace(': ', ':~')


class Section:

    def __init__(self, heading, options=None):
        self.heading = heading
        self.options = [] if options is None else options

    def __repr__(self):
        return '<{!s}: heading={!r}, len(options)={!r}>'.format(
            type(self).__name__, self.heading, len(self.options)
        )

    def __eq__(self, other):
        return all((
            self.heading == other.heading,
            self.options == other.options,
        ))


class Option:

    def __init__(self, name, argspec, desc):
        self.name = name
        self.argspec = argspec
        self.desc = desc

    def __repr__(self):
        return '<{!s}: name={!r}, argspec={!r}, desc={!r}>'.format(
            type(self).__name__, self.name, self.argspec, self.desc
        )

    def __eq__(self, other):
        return all((
            self.name == other.name,
            self.argspec == other.argspec,
            self.desc == other.desc,
        ))


def tokenise_boost_cpp_help(message):
    it = iter(message.split('\n'))
    next(it)  # skip a header line
    tokens = [
        ('empty', re.compile('^$')),
        ('section', re.compile('^(?P<name>.*):$')),
        ('option', re.compile(
            '^  (?P<name>\S+) (?P<argspec>( ?arg| ?\[.*?\]| ?\(.*?\))*)(?: +(?P<desc>.*))?$')),
        # intented text will be at least 3 spaces, exact number is
        # data-dependent
        ('indented', re.compile('^ {3,}(?P<text>.*)$')),
    ]
    for line in it:
        for name, regex in tokens:
            match = regex.match(line)
            if match:
                yield name, match
                break
        else:
            raise Exception('failed to match {!r}'.format(line))


def parse_boost_cpp_help(message):
    sections = []
    current_section = None
    current_option = None

    for name, match in tokenise_boost_cpp_help(message):
        if name == 'empty':
            continue
        elif name == 'section':
            current_section = Section(match.group('name'))
            sections.append(current_section)
            continue
        elif name == 'option':
            current_option = Option(**match.groupdict())
            current_section.options.append(current_option)
        elif name == 'indented':
            if current_option.desc is None:
                current_option.desc = match.group('text')
            else:
                current_option.desc = ' '.join(
                    (current_option.desc, match.group('text')))
        else:
            raise Exception(
                'unhandled token type {!r}: {!r}'.format(
                    name, match.string))

    return sections


def section_to_latex(sections):
    section_header = \
        '\\subsection{{{heading}}}\n' \
        '\\label{{subsection:{label}}}\n' \
        '\\begin{{tabular}}{{p{{0.3\\linewidth}}p{{0.2\\linewidth}}p{{0.45\\linewidth}}}}\n' \
        'Option & Arguments & Description \\\\ \\hline\n'
    option_line = '\\textbf{{{name}}} & \\textit{{{argspec}}} & {desc}\\\\\n'
    section_footer = '\\end{tabular}'
    return '\n\n'.join((
        section_header.format(
            heading=section.heading,
            label=section.heading.lower(),
        ) +
        ''.join((option_line.format(
            name=latex_encode(option.name).replace('-', '{-}'),
            argspec=latex_encode(option.argspec).replace('-', '{-}'),
            desc=latex_encode(option.desc),
        ) for option in section.options)) +
        section_footer
        for section in sections
    ))


def help_to_latex_main(argv):
    message = sys.stdin.read()
    sections = parse_boost_cpp_help(message)
    documentation = section_to_latex(sections)
    print(documentation)


if __name__ == '__main__':
    sys.exit(help_to_latex_main(sys.argv[:]))
