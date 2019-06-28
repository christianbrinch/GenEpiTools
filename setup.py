#!/usr/bin/env python
''' Setup file for GenEpiTools
'''

from setuptools import setup
import os

python_files = os.listdir('./genepitools/')
moduleNames = []
for file in python_files:
    if ".py" in file and not "init" in file and not ".swp" in file:
        moduleNames.append('genepitools.'+file.strip()[0:-3])

print(moduleNames)

setup(name='genepitools',
      version='0.2',
      description='Assorted Python tools for genomic epidemiology ',
      author='Christian Brinch',
      author_email='cbri@food.dtu.dk',
      py_modules=moduleNames
      )
