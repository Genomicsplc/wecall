# All content Copyright (C) 2018 Genomics plc
import logging
import subprocess


class ToolRunner(object):

    def __init__(self):
        self.return_code = None
        self.stdout = None
        self.stderr = None

    def log_output(self):
        log_output("stdout", self.stdout)
        log_output("stderr", self.stderr)
        logging.info('returncode: {}'.format(self.return_code))
        logging.info('')

    def run(self, command, cwd=None):
        # TODO: make this private
        proc = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd
        )

        self.stdout, self.stderr = proc.communicate()
        self.return_code = proc.returncode

        return self

    def start(self, command, cwd=None):
        log_command(command)
        self.run(command, cwd)
        self.log_output()
        return self


def subprocess_expectation_from_tool_runner(test_case, tool_runner):
    return SubprocessExpectation(
        test_case,
        tool_runner.stdout,
        tool_runner.stderr,
        tool_runner.return_code)


class SubprocessExpectation(object):
    def __init__(self, test_case, stdout, stderr, return_code):
        self.__test_case = test_case
        self.__stdout = stdout
        self.__stderr = stderr
        self.__return_code = return_code

    def failure(self):
        self.__test_case.assertNotEqual(0, self.__return_code)

    def success(self):
        self.__test_case.assertEqual(
            0, self.__return_code, msg="0!={}\nstderr:\n{!s}".format(
                self.__return_code, self.__stderr))


def log_command(command):
    logging.info("Running: `{}`".format(" ".join(command)))
    logging.info('')


def log_file(filename):
    if filename is None:
        return
    try:
        with open(filename) as fp:
            log_output(filename, fp.read())
    except IOError:
        logging.info('{!r} not found'.format(filename))


def log_output(title, output):
    logging.info("{}:".format(title))
    for line in str(output).split('\n'):
        logging.info('> ' + line)
    logging.info('')


def log_bam_file(reference_filename, bam_filename, chrom=None):
    command = ['samtools', 'tview', '-dT', bam_filename, reference_filename]
    if chrom is not None:
        print("Contig: {}".format(chrom))
        command.extend(["-p", chrom])
    log_output(bam_filename, subprocess.check_output(command))

    command_2 = ['samtools', 'view', bam_filename]
    if chrom is not None:
        print("Contig: {}".format(chrom))
    log_output(bam_filename, subprocess.check_output(command_2))
