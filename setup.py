#!/usr/bin/env python

import sys, os
from codecs import open

try:
    from setuptools import setup
except ImportError:   
    sys.stderr.write( "Could not import setuptools ! Please install it with: pip install setuptools " )
   
if sys.version_info < (2, 7):
    print "At least Python 2.7 is required! Please install Python 2.7 on your OS"
    exit(1)
   
here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='htseq-clip',
    version='0.1.14',
    description='htseq-clip: a pipeline for the analysis of iCLIP datasets',
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
	install_requires=['biopython', 'bokeh', 'HTSeq'],
	
    packages=['clip'],
	entry_points = {
        'console_scripts': ['htseq-clip=clip.command_line:main'],
    } 
)