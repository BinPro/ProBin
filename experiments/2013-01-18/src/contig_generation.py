#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser
from helpers import sample_contigs

from Bio import SeqIO

def main(genomes, n, min_length, max_length):
    for g in genomes:
        contigs = sample_contigs(g, n, min_length, max_length)
        sys.stdout.write("\n".join(contigs) + "\n\n")


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('files', nargs='*', 
        help='specify input files, default is stdin')
    parser.add_argument('-o', '--output', 
        help='specify the output file.  The default is stdout')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser.add_argument('-n', '--number', nargs='?', default=100, type=int,
        help='specify the number of contigs to generate from each genome, default is 100')
    parser.add_argument('--max_length', nargs='?', default=5000, type=int,
        help='specify the maximum contig length, default is 5000')
    parser.add_argument('--min_length', nargs='?', default=1000, type=int,
        help='specify the minimum contig length, default is 1000')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')

    handle = fileinput.input(args.files)
    genomes = list(SeqIO.parse(handle,"fasta"))
    handle.close()
    if args.verbose:
        sys.stderr.write("Number of genomes read: %i %s" % (len(genomes),os.linesep))
        sys.stderr.write("Number of contigs to be generated from each genome: %i %s" % (args.number, os.linesep) )
        sys.stderr.write("Max length of contigs to be generated from each genome: %i %s" % (args.max_length, os.linesep) )
        sys.stderr.write("Min length of contigs to be generated from each genome: %i %s" % (args.min_length, os.linesep) )
    
    main(genomes, args.number, args.min_length, args.max_length)
