#!/usr/bin/env python
import fileinput
import sys
import os
from argparse import ArgumentParser

from Bio import SeqIO

from probin.model.composition import  multinomial as ml

def main(contigs,kmer_len,verbose):
    signatures = ml.calculate_signatures(kmer_len, contigs)
    size_possible_kmers = 4**kmer_len 
    prob = [1.0/size_possible_kmers]*size_possible_kmers
    probability = ml.probability(signatures[0],prob)
    log_probability = ml.log_probability(signatures[0],prob)
    print probability
    print log_probability

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
    parser.add_argument('-m', '--model', nargs='?', default=1, type=int,
        help='specify the model number to use, default 1')
    args = parser.parse_args()
    
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')

    handle = fileinput.input(args.files)
    contigs = list(SeqIO.parse(handle,"fasta"))
    handle.close()
    if args.verbose:
        sys.stderr.write("Number of contigs read: %i %s" % (len(contigs),os.linesep))
    main(contigs,args.kmer, args.verbose)
