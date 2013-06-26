import mixture
from ProBin import _get_contigs
from probin.model.composition.log_signatures import signatures_to_log
import os
from os import getpid
from argparse import ArgumentParser
import numpy as np
from random import randint

def _initialize(contigs,K):
    rs = np.random.RandomState(seed=randint(0,10000)+getpid())
    N,D = contigs.shape
    ind = rs.choice(N,K,True)
    indx = np.arange(N)
    mu = np.zeros((K,D))
    for i,centroid in enumerate(ind):
        mu[i] = contigs[indx==centroid]
    
    sigmas = [np.identity(D)]*K
    return mu, sigmas

def main(K,contig_file,kmer_len,max_iter,epsilon):
    contigs,idx = _get_contigs(contig_file,kmer_len)

    log_contigs = signatures_to_log(contigs)
    N,D = log_contigs.shape
    mu,sigmas = _initialize(log_contigs,K)
    import sys
    sys.stderr.write(str(log_contigs) +'\n')
    comps = []
    for i in xrange(K):
        comps.append(mixture.MultiNormalDistribution(D,mu[i],sigmas[i]))

    m = mixture.MixtureModel(K,[1/float(K)]*K,comps)
    data = mixture.DataSet()
    data.fromArray(log_contigs)
    m.EM(data,max_iter,epsilon)

    return m

def main_parser():
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    parser.add_argument('composition_file',
                        help='Specify composition csv file')
    parser.add_argument('-o','--output',
                        help='Specify composition csv file')
    parser.add_argument('composition_file',
                        help='Specify composition csv file')
    
