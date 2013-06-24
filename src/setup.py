#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os

version = '0.2'

setup(name='probin',
      version=version,
      description="Probabilistic Binning",
      long_description="""\
To be done""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='Python Scilifelab Metagenomics Binning Clustering Contig',
      author='Brynjar Smari Bjarnason, Johannes Alneberg',
      author_email='brynjar.bjarnason@scilifelab.se',
      url='www.github.com/BinPro/ProBin',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      scripts=["ProBin.py"],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
