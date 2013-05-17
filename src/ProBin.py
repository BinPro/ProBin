#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script for clustering metagenomic contigs based on sequence composition and
correlation between many samples."""
import sys
import os

from Bio import SeqIO

from probin.dna import DNA
from probin.parser import main_parser
from probin.preprocess import main_preprocess

def main(contigs,model,clustering,cluster_count,verbose):
    uniform_prob = {}
    for i in xrange(DNA.kmer_hash_count):
        uniform_prob[i]= 1.0/float(DNA.kmer_hash_count)
    (clusters,clust_prob, centroids) = clustering.cluster(contigs, model.log_probability,model.fit_nonzero_parameters, cluster_count=cluster_count ,centroids=None, max_iter=100, repeat=10,epsilon=0.000001)
    return (clusters,clust_prob,centroids)


def write_clustering_result(clusters, cluster_evaluation, centroids, arguments, output):
    #CLUSTERING INFORMATION OUTPUT
    RESULT=["#{divide}",
            "#Clustering based on parameters: {args}",
            "#Result written to files starting with: {directory}",
            "#Clustering evaluation: {clust_prob}",
            "#<Cluster sizes>",
            "{cluster_freq}",
            "{clusters}\n"]
    repr_centroids = ["#Centroid {0},{1}".format(i,",".join(map(str,centroid))) for i,centroid in enumerate(centroids)]
    cluster_sizes = [len(c) for c in clusters]
    tot_c = float(sum(cluster_sizes))
    cluster_freq = ["#Cluster {0}:\t{1}\t{2}".format(i,c,c/tot_c) for i,c in enumerate(cluster_sizes)]
    cluster_contigs_id =  [ "Cluster {0}\t{1}".format(i,",".join([contig.id for contig in cluster]) )  for i,cluster in enumerate(clusters)]
    params =   {"args":arguments, "clust_prob":cluster_evaluation,
                "centroids":os.linesep.join(repr_centroids),
                "divide":"="*70,
                "directory":output,
                "cluster_freq":os.linesep.join(cluster_freq),
                "clusters":os.linesep.join(cluster_contigs_id)}
    with open(output,"w") as clustinf:
        clustinf.write(os.linesep.join(RESULT).format(**params))

def _get_contigs(arg_file):
    try:
        with open(arg_file) as handle:
            seqs = list(SeqIO.parse(handle,"fasta"))
    except IOError as error:    
        print >> sys.stderr, "Error reading file %s, message: %s" % (error.filename,error.message)
        sys.exit(-1)
    except Exception as error:
        print >> sys.stderr, "Error reading file %s, message: %s" % (error.filename,error.message)
        sys.exit(-1)

    contigs = [DNA(x.id, x.seq.tostring().upper(), calc_sign=True) for x in seqs]

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
            
        if args.verbose:
            print >> sys.stderr, "parameters: %s" % (args)
            print >> sys.stderr, "Reading file and generating contigs"

        if not os.path.isdir(os.path.abspath(args.output)):
            args.output = os.getcwd()
        output = os.sep.join([os.path.abspath(args.output),"{0}_k{1}_c{2}_{3}".format(os.path.basename(args.file),args.kmer,args.cluster_count,args.algorithm)])
        if os.path.isfile(output):
            from datetime import datetime
            output = "{0}_{1}".format(output,datetime.now().strftime("%Y-%m-%d-%H.%M"))
        print >> sys.stderr, "Result files created in: %s" % (os.path.dirname(output))
            

        if args.verbose:
            print >> sys.stderr, "parameters: %s" %(args)
        
        DNA.generate_kmer_hash(args.kmer)
    
        contigs = _get_contigs(args.file)
        
        (clusters,clust_prob,centroids) = main(contigs,model,algorithm,args.cluster_count, args.verbose)
    
        write_clustering_result(clusters,clust_prob,centroids,args,output)

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

