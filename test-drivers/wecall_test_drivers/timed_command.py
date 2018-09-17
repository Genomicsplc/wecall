# All content Copyright (C) 2018 Genomics plc
from wecall_test_drivers.tool_runner import ToolRunner
import json
import logging
import tempfile
import psutil
import time


class TimedCommand(ToolRunner):

    def __init__(self):
        ToolRunner.__init__(self)
        self.user_time = None
        self.system_time = None

    @property
    def times(self):
        return {"user_time": self.user_time, "system_time": self.system_time}

    def dump_timing_json(self, filename):
        with open(filename, "w") as json_fp:
            json.dump(self.times, json_fp, indent=4, sort_keys=True)
            json_fp.write("\n")

    def log_output(self):
        ToolRunner.log_output(self)
        logging.info("user_time: {}".format(self.user_time))
        logging.info("system_time: {}".format(self.system_time))

    def run(self, command, cwd=None):
        with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
            process = psutil.Popen(
                command, stdout=stdout, stderr=stderr, cwd=cwd)
            while process.status() != psutil.STATUS_ZOMBIE:
                time.sleep(0)
            stdout.seek(0)
            self.stdout = stdout.read()
            stderr.seek(0)
            self.stderr = stderr.read()

        times = process.cpu_times()
        self.user_time, self.system_time = times.user, times.system
        process.wait()
        self.return_code = process.returncode

        return self
