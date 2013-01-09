#!/usr/bin/python
# -*- coding: latin-1 -*-

"""
CRP will read files on fasta format and a dispersion parameter, and return a random binning of the contigs in the fasta file.

"""
import argparse
import sys
import Bio

def main():
  """
  Read fasta files and run Chinese Restaurant Process to create bins. Return the bins and elements

  """
  

if __name__=="__main__":
  parser = argparse.ArgumentParser(description='Read in a file and a dispersion parameter')
  parser.add_argument('infile', nargs='+', type=argparse.FileType('r'),default=sys.stdin)
  parser.add_argument('alpha', metavar='N', type=float, help='alpha, the dispersion parameter')
  args = parser.parse_args()
  main()
