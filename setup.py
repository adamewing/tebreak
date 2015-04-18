#!/usr/bin/env python

from setuptools import setup

setup(
    name='TEBreak',
    version='0.0.1',
    author='Adam Ewing',
    author_email='adamewing@gmail.com',
    description=("Insertion finder for high throughput sequence data"),
    license='MIT',
    url='https://github.com/adamewing/tebreak',
    packages=['tebreak'],
    install_requires = [
        'pysam>=0.8.1',
        'bx-python>=0.5.0',
        'scipy>=0.14.0',
        'numpy>=1.9.0',
        'align>=0.1',
    ],
    dependency_links = [
        'git+ssh://git@github.com:adamewing/align.git#egg=align-0.1'
    ],

)


