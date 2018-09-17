# All content Copyright (C) 2018 Genomics plc
import os
from wecall_test_drivers.wecall_config_file_test_runnner import WecallConfigFileRunnerTest


class TestsFromConfigFiles(WecallConfigFileRunnerTest):
    def test_basic_config(self):
        test_dir = os.path.dirname(__file__)
        fileName = os.path.join(test_dir, "config", "basic_calling_test.config")
        self.run_from_config_file(fileName)
