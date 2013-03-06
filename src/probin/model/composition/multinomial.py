#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from numpy import log

def fit_parameters(sig):
    par = {}
    n = sum(sig.values())
    for i,v in sig.items():
        par[i] = v/float(n)
    return par

def log_probability(signature, prob_vector):
    phi = sum(signature.values())
    log_prod = 0
    for i,cnt in signature.iteritems():
        denom = _log_fac(cnt)
        log_prod += (log(prob_vector[i])*cnt) - denom
    return log_prod + _log_fac(phi)

def _log_fac(i):
    return sum(log(range(1,i+1)))
