#!/usr/bin/env python

import sys
import subprocess
from setuptools import setup, find_packages

def check_bwa():
    p = subprocess.Popen(['bwa'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()
        if line.startswith('Version:'):
            major, minor, sub = line.strip().split()[1].split('.')
            sub = sub.split('-')[0]
            if int(major) >= 0 and int(minor) >= 7 and int(sub) >= 12:
                return True
    return False


def check_samtools():
    p = subprocess.Popen(['samtools'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()
        if line.startswith('Version:'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False


def check_minia():
    p = subprocess.Popen(['minia'], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode()
        if line.startswith('[minia options]'):
            return True
    return False


def check_LAST():
    p = subprocess.Popen(['lastal'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()
        if line.startswith('lastal'):
            return True 
    return False


def check_python():
    return sys.hexversion >= 0x20702f0    


def check_exonerate():
    p = subprocess.Popen(['exonerate'], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode()
        if line.startswith('exonerate from exonerate'):
            major, minor = line.strip().split()[-1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 2 and int(minor) >= 2:
                return True
    return False


if __name__ == '__main__':
    if not check_python(): sys.exit('Dependency problem: python >= 2.7.2 is required')
    if not check_bwa(): sys.exit('Dependency problem: bwa >= 0.7.12 not found')
    if not check_samtools(): sys.exit('Dependency problem: samtools >= 1.0 not found')
    if not check_minia(): sys.exit('Dependency problem: minia not found')
    if not check_LAST(): sys.exit('Dependency problem: LAST >= 548 not found')
    if not check_exonerate(): sys.exit('Dependency problem: exonerate >= 2.2 not found (required for post-filter)')

setup(
    name='TEBreak',
    version='1.1',
    author='Adam Ewing',
    author_email='adam.ewing@gmail.com',
    description=("Insertion finder for high throughput sequence data"),
    license='MIT',
    url='https://github.com/adamewing/tebreak',
    scripts=['tebreak/tebreak'],
    packages=find_packages(),
    install_requires = [
        'pysam>=0.8.1',
        'bx-python>=0.5.0',
        'scipy>=0.14.0',
        'numpy>=1.9.0',
        'scikit-bio>=0.5.5',
    ]

)


