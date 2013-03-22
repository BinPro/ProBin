#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script for clustering metagenomic contigs based on sequence composition and
correlation between many samples."""
import fileinput
import sys
import os
from argparse import ArgumentParser

from Bio import SeqIO

from probin.dna import DNA

def main(contigs,model,clustering,verbose):
    uniform_prob = {}
    for i in xrange(DNA.kmer_hash_count):
        uniform_prob[i]= 1.0/float(DNA.kmer_hash_count)
    (clust_prob, centroids, clusters) = clustering.cluster(contigs, model, cluster_count=3 ,centroids=None, max_iter=100, repeat=10)
    
    return (clust_prob,centroids,clusters)


def print_clustering_result(cluster_evaluation, centroids, clusters, arguments):
    RESULT="""#Clustering based on parameters: {args}.
#clustering evaluation: {clust_prob}
#<Centroids>
#{centroids}
{clusters}"""
    c = [">Cluster {0}{1}{2}".format(i,os.linesep, os.linesep.join(map(str,cluster))) for i,cluster in enumerate(clusters)]
    params =   {"args":arguments, "clust_prob":cluster_evaluation,\
                "clusters":os.linesep.join(c),
                "centroids":" ".join(map(str,centroids[i]))}
    print RESULT.format(**params)

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
    parser.add_argument('-c', '--clustering', default='kmeans', type=str, choices=['kmeans'],
        help='specify the clustering to use, default kmeans.')
    parser.add_argument('-cc', '--cluster_count', default=10, type=int,
        help='specify the number of cluster to use')
    args = parser.parse_args()
    
    try:
        model = __import__("probin.model.composition.{0}".format(args.model_composition),globals(),locals(),["*"],-1)
    except ImportError:
        print "Failed to load module {0}. Will now exit".format(args.model_composition)
        sys.exit(-1)
        
    try:
        clustering = __import__("probin.binning.{0}".format(args.clustering),globals(),locals(),["*"],-1)
    except ImportError:
        print "Failed to load module {0}. Will now exit".format(args.clustering)
        sys.exit(-1)
        
    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
    
    if args.verbose:
        print >> sys.stderr, "parameters: %s" % (args)
        print >> sys.stderr, "Reading file and generating contigs"
    DNA.generate_kmer_hash(args.kmer)
    try:
        handle = fileinput.input(args.files)
        contigs = [DNA(x.id, x.seq.tostring().upper(),calc_sign=True) for x in list(SeqIO.parse(handle,"fasta"))]
    except IOError as error:
        print >> sys.stderr, "Error reading file %s, message: %s" % (error.filename,error.message)
        sys.exit(-1)
    finally:
        handle.close()

    if args.verbose:
        print >> sys.stderr, "parameters: %s" %(args)
    
    (clust_prob,centroids,clusters) = main(contigs,model,clustering, args.verbose)
    print_clustering_result(clust_prob,centroids,clusters,args)