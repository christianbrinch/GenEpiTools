#!/usr/bin/env python
''' Setup file for GenEpiTools
'''

from distutils.core import setup

setup(name='genepitools',
      version='0.01',
      description='Assorted Python tools for genomic epidemiology ',
      author='Christian Brinch',
      author_email='cbri@food.dtu.dk',
      py_modules=['__init__', 'metadata', 'functions'])
