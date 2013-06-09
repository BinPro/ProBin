#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from scipy.special import gammaln
import numpy as np
from collections import Counter
import sys
def fit_parameters(dna_l):
    sig = Counter()
    [sig.update(part.signature) for part in dna_l]
    par = np.zeros(dna_l[0].kmer_hash_count)
    for key,cnt in sig.iteritems():
        par[key] += cnt
    par /= np.sum(par)
    return par
def fit_nonzero_parameters(composition,expected_clustering=None):
    if expected_clustering == None:
        expected_clustering = np.ones((len(composition),1))
    pseudo_sig = (expected_clustering.T).dot(composition) + 1
    pseudo_sig /= np.sum(pseudo_sig,axis=1,keepdims=True)
    if len(pseudo_sig) == 1:
        pseudo_sig = pseudo_sig[0]
    return pseudo_sig
def log_probability(seq, prob_vector):
    signature_vector = np.zeros(np.shape(prob_vector))

    for key,value in seq.signature.iteritems():
        signature_vector[key] = value

    return np.sum((signature_vector * np.log(prob_vector)) - _log_fac(signature_vector)) + _log_fac(np.sum(signature_vector))
def log_probabilities(composition, prob_vectors):
    if len(prob_vectors.shape) == 1:
        prob_vectors = np.array([prob_vectors])    
    log_qs = np.zeros((len(composition),len(prob_vectors)))
    for i,sign in enumerate(composition):
        log_qs[i] = np.sum((sign * np.log(prob_vectors)) - _log_fac(sign),axis=1) + _log_fac(np.sum(sign))
    return log_qs
def _log_fac(i):
    # gammaln produces the natural logarithm of the factorial of i-1
    return gammaln(i+1)
