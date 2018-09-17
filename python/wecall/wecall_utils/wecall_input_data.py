# All content Copyright (C) 2018 Genomics plc
class InputData(object):

    def __init__(self, tags, filenames):
        self.tags = tags
        self.__filenames = filenames

    @property
    def filenames(self):
        return self.__filenames


class WecallInputData(InputData):

    def __init__(self, bam_filenames, reference_filename):
        InputData.__init__(self, set(), set())
        self.bam_filenames = bam_filenames
        self.reference_filename = reference_filename

    @property
    def filenames(self):
        filenames = set()
        for bam_filename in self.bam_filenames:
            filenames.update({bam_filename, bam_filename + ".bai"})
        filenames.update(
            {self.reference_filename, self.reference_filename + ".fai"})
        return filenames
