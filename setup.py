# -*- encoding: utf-8 -*-

import io
import os
import re
from setuptools import setup
import sys

if sys.version_info[0] < 2 or \
        sys.version_info[0] == 2 and sys.version_info[1] < 7:
    sys.stderr.write('Error in seqpoet setup\n')
    sys.stderr.write('You need at least version 2.7 of Python to use seqpoet\n')
    sys.exit(1)

if sys.version_info[0] >= 3:
    sys.stderr.write('Error in seqpoet setup\n')
    sys.stderr.write('This package only works with Python 2 at the moment\n')
    sys.stderr.write('Please use Python 2.x, x >= 7\n')
    sys.exit(1)

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(name='seqpoet',
    version=find_version('seqpoet/__init__.py'),
    description='In silico PCR and operon extraction',
    url='https://github.com/maehler/seqpoet',
    author='Niklas Mähler',
    author_email='niklas.mahler@gmail.com',
    maintainer='Niklas Mähler',
    maintainer_email='niklas.mahler@gmail.com',
    license='MIT',
    packages=['seqpoet'],
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'],
    scripts=['bin/seqpoet'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2 :: Only',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
