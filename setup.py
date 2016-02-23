"""
__author__ = "Marko Fritz"
__copyright__ = "Copyright 2016, Marko Fritz"
__email__ = "marko.fritz@embl.de"
__license__ = "EMBL"
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
    version='0.1.0',
    description='htseqCLIP - a pipeline for the analysis of iCLIP datasets',
	long_description=long_description,
    url='https://bitbucket.org/htseq-clip/htseq-clip',
    author='Marko Fritz',
    author_email='marko.fritz@embl.de',
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
    packages=['htseq-clip', 'htseq-clip/lib'],
    install_requires=['biopython', 'bokeh', 'HTSeq'],
    entry_points={
        'console_scripts': [
            'htseq-clip=htseqclip:main',
        ],
    },
)