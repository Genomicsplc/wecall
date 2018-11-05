# All content Copyright (C) 2018 Genomics plc
import os
import subprocess
from wecall.common.exceptions import weCallException
from wecall.genomics.reference_genome import InMemoryReferenceGenome
from wecall_test_drivers.tool_runner import ToolRunner


class FastaFileBuilder(object):

    def __init__(self, filename, line_length=80):
        assert(line_length > 0)
        self.filename = filename
        self.line_length = line_length
        self.__reference_genome = InMemoryReferenceGenome()

    def reference_genome(self):
        return self.__reference_genome

    def with_chrom(self, chrom, sequence, pos_from=0):
        return self.__reference_genome.with_chrom(chrom, sequence, pos_from)

    def build(self):
        with open(self.filename, "w") as fasta_file:
            for chrom_name in self.__reference_genome.chromosomes():
                sequence = self.__reference_genome.fetch(chrom_name)
                fasta_file.write(">{}\n".format(chrom_name))
                for offset in range(0, len(sequence), self.line_length):
                    line = sequence[offset:offset + self.line_length] + '\n'
                    fasta_file.write(line)

        return self

    def index(self):
        tool_runner = ToolRunner()
        tool_runner.start(
            [os.path.join(os.environ['WECALL_BIN'], "samtools"), "faidx", self.filename])

        if tool_runner.return_code != 0:
            raise weCallException("")
        else:
            return self
