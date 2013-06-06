#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Script for clustering metagenomic contigs based on sequence composition and
correlation between many samples."""
import sys
import pandas as p # Used by _get_coverage

from Bio import SeqIO

from probin.dna import DNA
from probin.output import Output
from probin.parser import main_parser
from probin.preprocess import main_preprocess


def main(contigs,model_composition,algorithm,cluster_count,verbose,iterations,repeat,epsilon,**kwargs):
    if kwargs['coverage'] is None:
        (clusters,clust_prob, centroids) = algorithm.cluster(contigs, model_composition.log_probabilities ,model_composition.fit_nonzero_parameters, cluster_count=cluster_count ,centroids=None, max_iter=iterations, repeat=repeat,epsilon=epsilon, verbose=verbose, **kwargs)
    else:
        (clusters,clust_prob, centroids) = algorithm.cluster(contigs, kwargs['model_coverage'], cluster_count=cluster_count ,centroids=None, max_iter=iterations, repeat=repeat,epsilon=epsilon,verbose=verbose,**kwargs)
    return (clusters,clust_prob,centroids)

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

def _get_coverage(arg_file):
    try:
        return p.io.parsers.read_table(arg_file,sep='\t',index_col=0)
    except Exception as error:
        print >> sys.stderr, "Error reading file %s, message: %s" % (error.filename,error.message)
        sys.exit(-1)

if __name__=="__main__":
    parser = main_parser()
    args = parser.parse_args()
    
    if args.script == 'probin':
        #Import composition model and prep data for that
        if args.model_composition is not None:
            try:
                model_composition = __import__("probin.model.composition.{0}".format(args.model_composition),globals(),locals(),["*"],-1)
                DNA.generate_kmer_hash(args.kmer)
                contigs = _get_contigs(args.file)

            except ImportError:
                print "Failed to load module {0}. Will now exit".format(args.model_composition)
                sys.exit(-1)
        else:
            model_composition=None
        
        #Import coverage model and prep data for that
        if args.model_coverage is not None:
            try:
                model_coverage = __import__("probin.model.coverage.{0}".format(args.model_coverage),globals(),locals(),["*"],-1)
                coverage = _get_coverage(args.coverage_file)

            except ImportError:
                print "Failed to load module {0}. Will now exit".format(args.model_coverage)
                sys.exit(-1)
        else:
            model_coverage = None
            coverage = None

        #Import clustering algorithm
        try:
            algorithm = __import__("probin.binning.{0}".format(args.algorithm),globals(),locals(),["*"],-1)
        except ImportError:
            print "Failed to load module {0}. Will now exit".format(args.algorithm)
            sys.exit(-1)

        #Prep output settings
        Output.set_output_path(args.output,args.file,args.kmer,args.cluster_count,args.algorithm)
        
    


        (clusters,clust_prob,centroids) = main(contigs,model_composition,algorithm,args.cluster_count, args.verbose, args.iterations, args.repeat,args.epsilon, model_coverage=model_coverage,coverage=coverage, last_data=args.last_data,first_data=args.first_data,read_length=args.read_length)

        Output.write_clustering_result(clusters,clust_prob,centroids,args)

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

