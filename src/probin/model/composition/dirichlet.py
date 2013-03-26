"""Implementation of a dirichlet model based on sequence composition"""

from numpy import log
import numpy as np
from scipy.special import gammaln
from scipy.optimize import fmin_l_bfgs_b

def fit_nonzero_parameters(dna_l,kmer_hash_count):
    alpha0 = [1.0]*kmer_hash_count
    alpha_bounds = [(0.0,None)]*kmer_hash_count
    alpha_fit = fmin_l_bfgs_b(neg_log_probability_l,alpha0,args=(dna_l,),bounds=alpha_bounds, approx_grad=True)
    return np.array(alpha_fit[0])

def log_probability(seq,alpha):
    pc = np.fromiter(seq.pseudo_counts,np.dtype('u4'),count=alpha.shape[0])
    return _log_probability_body(pc,alpha)

def log_probability_test(pseudo_counts,alpha):
    # pseudo_counts given directly mainly for testing purposes
    pc = np.array(pseudo_counts)
    return _log_probability_body(pc,alpha)

def _log_probability_body(pc,alpha):
    A = np.sum(alpha)
    N = np.sum(pc)
    s = np.sum(gammaln(pc+alpha)-gammaln(alpha))
    return gammaln(A) - gammaln(A+N) + s 


def neg_log_probability_l(alpha,dna_l):
    return -sum([log_probability(dna,alpha) for dna in dna_l])

