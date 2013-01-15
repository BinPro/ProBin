#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser

from Bio import SeqIO


def main(contigs,kmer_len,verbose):
    print kmer_len


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
        help='specify the output file.  The default is stdout')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser.add_argument('-k', '--kmer', nargs='?', default=4, type=int,
        help='specify the length of kmer to use, default 4')
    args = parser.parse_args()
    
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')

    handle = fileinput.input(args.files)
    contigs = list(SeqIO.parse(handle,"fasta"))
    handle.close()
    if args.verbose:
        sys.stderr.write("Number of contigs read: %i %s" % (len(contigs),os.linesep))
    main(contigs,args.kmer, args.verbose)
