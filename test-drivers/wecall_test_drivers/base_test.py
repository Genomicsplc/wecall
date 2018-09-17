# All content Copyright (C) 2018 Genomics plc
from logging import FileHandler, StreamHandler, DEBUG, INFO, getLogger, Formatter
import shutil
import unittest
import os
import sys


class BaseTest(unittest.TestCase):

    def setUp(self):
        self.work_dir = os.path.join(
            os.environ["ECHIDNA_TEST_RESULTS"],
            *self.id().split("."))
        if os.path.exists(self.work_dir):
            shutil.rmtree(self.work_dir)
        os.makedirs(self.work_dir)

        logger = getLogger()
        logger.setLevel(DEBUG)
        fh = FileHandler(os.path.join(self.work_dir, "test.log"))
        ch = StreamHandler(sys.stdout)
        logger.addHandler(configure_log_handler(fh, DEBUG))
        logger.addHandler(configure_log_handler(ch, INFO))

    def tearDown(self):
        logger = getLogger()
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)


def configure_log_handler(handler, level):
    handler.setLevel(level)
    formatter = Formatter('%(message)s')
    handler.setFormatter(formatter)
    return handler
