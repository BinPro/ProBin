#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from numpy import log

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
    phi = sum(signature.values())
    log_prod = 0
    for i,cnt in signature.items():
        denom = _log_fac(cnt)
        log_prod += (log(prob_vector[i])*cnt) - denom
    return log_prod + _log_fac(phi)

def _log_fac(i):
    return sum(log(range(1,i+1)))
