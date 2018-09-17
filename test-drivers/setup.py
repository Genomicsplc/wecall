# All content Copyright (C) 2018 Genomics plc
# -*- coding: utf8 -*-

import setuptools

setuptools.setup(
    name="wecall-test-drivers",
    url="www.genomicsplc.com",
    author="Genomics",
    author_email="help@genomicsplc.com",
    description="wecall-test-drivers",
    license="Genomics PLC Proprietary License",
    keywords="wecall-test-drivers",
    packages=setuptools.find_packages(),
    py_modules=['wecall-test-drivers'],
    install_requires=[(
        'pytest', 'pytest-cov', 'pytest-flakes', 'pytest-pep8', 'pytest-xdist',
        'testfixtures', 'pysam', 'psutil')],
)
