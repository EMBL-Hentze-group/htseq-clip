"""
__author__ = "Marko Fritz, Thomas Schwarzl"
__copyright__ = "Copyright 2016, Marko Fritz, Thomas Schwarzl"
__email__ = "marko.fritz@embl.de, schwarzl@embl.de"
__license__ = "MIT"
"""

from setuptools import setup
from codecs import open
from os import path
import sys

if sys.version_info < (2, 7):
    print "At least Python 2.7 is required!"
    exit(1)

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='htseq-clip',
    version='0.1.6',
    description='HTSeq clip - a pipeline for the analysis of iCLIP datasets',
	long_description=long_description,
    url='https://bitbucket.org/htseq-clip/htseq-clip',
    author='Marko Fritz, Thomas Schwarzl',
    author_email='marko.fritz@embl.de, schwarzl@embl.de',
	zip_safe=False,
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Environment :: Console',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ],
    packages=['clip'],
	scripts = ['scripts/htseq-clip'],
    install_requires=['biopython', 'bokeh', 'HTSeq'],
)