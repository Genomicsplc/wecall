# All content Copyright (C) 2018 Genomics plc
# -*- coding: utf8 -*-

import setuptools

setuptools.setup(
    name="wecall",
    url="www.genomicsplc.com",
    author="Genomics",
    author_email="help@genomicsplc.com",
    description="wecall",
    license="Genomics PLC Proprietary License",
    keywords="wecall",
    packages=setuptools.find_packages(),
    py_modules=['wecall'],
    install_requires=['pysam']
)
