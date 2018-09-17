# All content Copyright (C) 2018 Genomics plc
import os
import subprocess

ECHIDNA_BIN = os.environ["ECHIDNA_BIN"]


class TabixIndexer(object):

    def __init__(self, filename, file_type=None):
        self.filename = filename
        self.file_type = file_type

    @property
    def compressed_filename(self):
        return self.filename + ".gz"

    @property
    def compressed_filename_index(self):
        return self.compressed_filename + ".tbi"

    def bgzip(self):
        subprocess.call(
            [os.path.join(ECHIDNA_BIN, "bgzip"), "-f", self.filename])
        return self

    def index(self):
        self.bgzip()
        tabix_args = [os.path.join(ECHIDNA_BIN, "tabix"), "-f", ]
        if self.file_type == "VARINFO":
            tabix_args += ['-s', '1', '-b', '2', '-e', '3']
        elif self.file_type is not None:
            tabix_args += ["-p", self.file_type]
        tabix_args.append(self.compressed_filename)
        subprocess.check_call(tabix_args)
        return self
