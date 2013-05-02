#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from scipy.special import gammaln
import numpy as np
from collections import Counter

def fit_parameters(dna_l):
    sig = Counter()
    [sig.update(part.signature) for part in dna_l]
    par = np.zeros(dna_l[0].kmer_hash_count)
    for key,cnt in sig.iteritems():
        par[key] += cnt
    par /= np.sum(par)
    return par

def fit_nonzero_parameters(dna_l,expected_clustering=None):
    pseudo_sig = np.zeros((len(dna_l),dna_l[0].kmer_hash_count))
    for i,dna in enumerate(dna_l):
        pseudo_sig[i,:] = np.fromiter(dna.pseudo_counts,dtype=np.int) - 1 
    if expected_clustering == None:
        expected_clustering = np.ones((1,len(dna_l)))
    pseudo_sig = expected_clustering.dot(pseudo_sig).flatten()
    pseudo_sig += 1
    pseudo_sig /= np.sum(pseudo_sig)
    return pseudo_sig

def log_probability(seq, prob_vector):
    signature_vector = np.zeros(np.shape(prob_vector))

    for key,value in seq.signature.iteritems():
        signature_vector[key] = value

    return np.sum((signature_vector * np.log(prob_vector)) - _log_fac(signature_vector)) + _log_fac(np.sum(signature_vector))

def _log_fac(i):
    # gammaln produces the natural logarithm of the factorial of i-1
    return gammaln(i+1)
