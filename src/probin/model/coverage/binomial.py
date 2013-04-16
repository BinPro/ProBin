"""Implementation of a binomial model based on abundance variation."""

from numpy import log
import numpy as np
from scipy.special import gammaln

def fit_nonzero_parameters(dna_l):
    pass

def neg_log_probability_l(alpha,dna_l):
    pass

def log_probability(seq,q,M,L):
    """ Calculates the log probability for a single sequence
    
    seq - dna instance
    q   - vector with relative frequency of 
    frequency of cluster k in sample l
    """
    # q_k_l is the frequency of cluster k in sample l
    q_prime = q*len(seq)/float(L)
    M_i_l
    

def probability(seq,q,M, L):
    Q = np.prod(
