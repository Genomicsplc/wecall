# All content Copyright (C) 2018 Genomics plc
import os
import unittest

import pep8


class TestPep8(unittest.TestCase):
    """Run PEP8 on all files in this directory and subdirectories."""

    def setUp(self):
        base_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))))
        self.base_dirs = [os.path.join(base_dir, sub_dir) for sub_dir
                          in ["python", "scripts", "test", "test_drivers"]]

    def test_pep8(self):
        style = pep8.StyleGuide()
        style.options.max_line_length = 120  # because it isn't 1928 anymore
        errors = 0

        for base_dir in self.base_dirs:
            for root, _, files in os.walk(base_dir):
                python_files = [f for f in files if f.endswith('.py')]
                for pf in python_files:
                    check = style.check_files([os.path.join(root, pf)])
                    errors += check.file_errors
        self.assertEqual(errors, 0)
