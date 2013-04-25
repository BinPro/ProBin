#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script for clustering metagenomic contigs based on sequence composition and
correlation between many samples."""
import fileinput
import sys
import os
from itertools import izip
from argparse import ArgumentParser

from Bio import SeqIO

from probin.dna import DNA
from probin.binning import statistics as stats
from probin.parser import main_parser
from probin.preprocess import main_preprocess

def main(contigs,model,clustering,cluster_count,verbose):
    uniform_prob = {}
    for i in xrange(DNA.kmer_hash_count):
        uniform_prob[i]= 1.0/float(DNA.kmer_hash_count)
    (clusters,clust_prob, centroids) = clustering.cluster(contigs, model, cluster_count=cluster_count ,centroids=None, max_iter=100, repeat=10)
    return (clusters,clust_prob,centroids)


def print_clustering_result(clusters, cluster_evaluation, centroids, arguments, output):
    RESULT=["#{divide}",
            "#Clustering based on parameters: {args}",
            "#clustering evaluation: {clust_prob}",
            "#<Centroids>",
            "{centroids}",
            "{clusters}"]
    repr_centroids = ["#Centroid {0},{1}".format(i,",".join(map(str,centroid))) for i,centroid in enumerate(centroids)]
    
    c = [">Cluster {0}{1}{2}".format(i,os.linesep, os.linesep.join(map(str,cluster))) for i,cluster in enumerate(clusters)]
    params =   {"args":arguments, "clust_prob":cluster_evaluation,
                "clusters":os.linesep.join(c),
                "centroids":os.linesep.join(repr_centroids),
                "divide":"="*20}
    print>>output, os.linesep.join(RESULT).format(**params)

def _get_contigs(arg_files):
    try:
        handle = fileinput.input(arg_files)
        seqs = list(SeqIO.parse(handle,"fasta"))
    except IOError as error:    
        print >> sys.stderr, "Error reading file %s, message: %s" % (error.filename,error.message)
        sys.exit(-1)
    finally:
        handle.close()

    contigs = [DNA(x.id, x.seq.tostring().upper(), calc_sign=True) for x in seqs]
    try:
        for contig,seq in izip(contigs,seqs):
            contig.phylo = seq.description.split(" ",1)[1]
    except Exception as error:
        print >> sys.stderr, "No phylo information %s, message: %s" % (error,error.message)

    return contigs

if __name__=="__main__":
    parser = main_parser()
    args = parser.parse_args()
    
    if args.script == 'probin':
        try:
            model = __import__("probin.model.composition.{0}".format(args.model_composition),globals(),locals(),["*"],-1)
        except ImportError:
            print "Failed to load module {0}. Will now exit".format(args.model_composition)
            sys.exit(-1)
        try:
            algorithm = __import__("probin.binning.{0}".format(args.algorithm),globals(),locals(),["*"],-1)
        except ImportError:
            print "Failed to load module {0}. Will now exit".format(args.algorithm)
            sys.exit(-1)
        

        if args.output and args.output != '-':
            output = open(args.output, 'w')
        else:
            output = sys.stdout
    
        if args.verbose:
            print >> sys.stderr, "parameters: %s" % (args)
            print >> sys.stderr, "Reading file and generating contigs"
        
        DNA.generate_kmer_hash(args.kmer)
    
        contigs = _get_contigs(args.files)

        if args.verbose:
            print >> sys.stderr, "parameters: %s" %(args)

        (clusters,clust_prob,centroids) = main(contigs,model,algorithm,args.cluster_count, args.verbose)
        try:
            stats.get_statistics(contigs,clusters,output)
        except AttributeError:
            print >> sys.stderr, "Phylo attribute not available"
    
        print_clustering_result(clusters,clust_prob,centroids,args,output)
        output.close()

    elif args.script == 'preprocess':
        if args.output:
            args.output = open(args.output,'w+')
        else:
            args.output = sys.stdout
        if not args.contigs:
            print "Contigs file not correctly supplied, will now exit"
            sys.exit(-1)
        main_preprocess(args)
    else:
        pass

