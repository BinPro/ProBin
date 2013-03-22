#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from scipy.special import gammaln
import numpy as np

def fit_parameters(sig):
    par = {}
    n = sum(sig.values())
    for i,v in sig.items():
        par[i] = v/float(n)
    return par

def fit_nonzero_parameters(sig,kmer_hash_count):
    pseudo_sig = np.ones(kmer_hash_count)
    for key,cnt in sig.iteritems():
        pseudo_sig[key] += cnt
    pseudo_sig /= np.sum(pseudo_sig)
    return pseudo_sig

def log_probability(signature, prob_vector):
    signature_vector = np.zeros(np.shape(prob_vector))

    for key,value in signature.iteritems():
        signature_vector[key] = value

    return np.sum((signature_vector * np.log(prob_vector)) - _log_fac(signature_vector)) + _log_fac(np.sum(signature_vector))

def _log_fac(i):
    # gammaln produces the natural logarithm of the factorial of i-1
    return gammaln(i+1)
