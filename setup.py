#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
from setuptools import setup, find_packages

# Get the long description from the relevant file
with open('README.md') as f:
    readme = f.read()

setup(
    name='MicroscPSF-Py',
    version=0.1,
    description='Gibson-Lanni PSF calculation code.',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Kyle Douglass, Hazen Babcock',
    author_email='hbabcock@mac.com',
    url='https://github.com/MicroscPSF/MicroscPSF-Py',

    zip_safe=False,
    packages=find_packages(),

    package_data={},
    exclude_package_data={},
    include_package_data=True,

    requires=[],

    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

    license='MIT',
    keywords='PSF,microscopy',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        "Programming Language :: C",
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',        
    ],
)
