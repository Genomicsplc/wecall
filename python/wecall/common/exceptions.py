# All content Copyright (C) 2018 Genomics plc
"""
Basic exception classes to be used throughout the weCall code.
"""


class weCallException(Exception):
    """
    Base class for all exceptions. Everything we throw
    should derive from this.
    """

    def __init__(self, value):
        self.value = value
        self.message = value

    def __str__(self):
        return repr(self.value)


class weCallRuntimeException(weCallException):

    def __init__(self, return_code, result):
        self.return_code = return_code
        self.result = result
        self.value = "weCall exited with non-zero exit code {} for `{}`.".format(
            self.return_code, self.result)
