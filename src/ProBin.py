#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script for clustering metagenomic contigs based on sequence composition and
correlation between many samples."""
import fileinput
import sys
import os
from argparse import ArgumentParser

from Bio import SeqIO

from probin.model.composition import  multinomial as ml
from probin.dna import DNA

def main(contigs,verbose):
    size_possible_kmers = 4**DNA.kmer_len
    uniform_prob = {}
    for i in xrange(DNA.kmer_hash_count):
        uniform_prob[i]= 1.0/float(DNA.kmer_hash_count)
    for contig in contigs:
        log_probability = ml.log_probability(contig.signature,uniform_prob)
        print log_probability



if __name__=="__main__":
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    parser.add_argument('files', nargs='*', 
        help='specify input files on FASTA format, default is stdin')
    parser.add_argument('-o', '--output', 
        help='specify the output file.  The default is stdout')
    parser.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser.add_argument('-k', '--kmer', default=4, type=int,
        help='specify the length of kmer to use, default 4')
    parser.add_argument('-mc', '--model_composition', default='multinomial', type=str, choices=['multinomial'],
        help='specify the composition model to use, default multinomial.')
    args = parser.parse_args()
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
    
    if args.verbose:
        print >> sys.stderr, "parameters: %s" % (args)
        print >> sys.stderr, "Reading file and generating contigs"
    DNA.generate_kmer_hash(args.kmer)
    try:
        handle = fileinput.input(args.files)
        contigs = [DNA(x.id, x.seq.tostring().upper()) for x in list(SeqIO.parse(handle,"fasta"))]
    except IOError as error:
        print >> sys.stderr, "Error reading file %s, message: %s" % (error.filename,error.message)
    finally:
        handle.close()

    if args.verbose:
        print >> sys.stderr, "parameters: %s" %(args)
    main(contigs, args.verbose)
