"""Implementation of a dirichlet model based on sequence composition"""

from numpy import log
import numpy as np
from scipy.special import gammaln
from scipy.optimize import fmin_l_bfgs_b

def fit_nonzero_parameters(dna_l,kmer_hash_count):
    alpha0 = [1.0]*kmer_hash_count
    alpha_bounds = [(0.0,None)]*kmer_hash_count
    alpha_fit = fmin_l_bfgs_b(neg_log_probability_l,alpha0,args=(dna_l,),bounds=alpha_bounds, approx_grad=True)
    return alpha_fit[0]

def log_probability(seq,alpha):
    A = sum(alpha)
    N = sum(seq.pseudo_counts)
    D = len(alpha)
    return gammaln(A) - gammaln(A+N) \
        + sum([gammaln(seq.pseudo_count(j)+alpha[j]) \
                   - gammaln(alpha[j]) for j in xrange(D)])

def neg_log_probability_l(alpha,dna_l):
    return -sum([log_probability(dna,alpha) for dna in dna_l])

