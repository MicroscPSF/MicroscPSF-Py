#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from setuptools import setup, find_packages

version = "0.0"
description = "Gibson-Lanni PSF calculation code."
long_description = ""

setup(
    name='storm_analysis',
    version=version,
    description=description,
    long_description=long_description,
    author='?',
    author_email='?',
    url='https://github.com/?',

    zip_safe=False,
    packages=find_packages(),

    package_data={},
    exclude_package_data={},
    include_package_data=True,

    requires=[],

    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    
    license="",  
    keywords='PSF,microscopy',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: ?',
        "Programming Language :: C",
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',        
    ],
)
