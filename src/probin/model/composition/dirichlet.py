"""Implementation of a dirichlet model based on sequence composition"""

from numpy import log
import numpy as np
from scipy.special import gammaln

def fit_nonzero_parameters(sig_l,kmer_hash_count):
    pass

def log_probability(sig,alpha):
    return log_density(sig,alpha)

def log_density(sig, alpha):
    A = sum(alpha)
    N = sum(sig)
    D = len(sig)
    return gammaln(A) - gammaln(A+N) \
        + sum([gammaln(sig[j]+alpha[j]) - gammaln(alpha[j]) for \
                   j in xrange(D)])

def log_density_l(sig_l,alpha):
    return sum([log_probability(sig,alpha) for sig in sig_l])

