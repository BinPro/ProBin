#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from numpy import log
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
    log_prod = _log_fac(sum(signature.values()))
    log_prob_vector = log(prob_vector)
    
    signature_vector = np.zeros(log_prob_vector.shape)
    
    for key,value in signature.iteritems():
        signature_vector[key] = value
    log_prod = sum(signature_vector*log_prob_vector - _log_fac(signature_vector))

    for i,cnt in signature.iteritems():
        denom = _log_fac(cnt)
        log_prod += (log_prob_vector[i]*cnt) - denom
    return log_prod

def _log_fac(i):
    # gammaln produces the natural logarithm of the factorial of i-1
    return gammaln(i+1)
