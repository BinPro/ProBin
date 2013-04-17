"""Implementation of a binomial model based on abundance variation."""

from numpy import log
import numpy as np
from scipy.special import gammaln

def fit_nonzero_parameters(dna_l,M):
    """ 
    Estimates pseudo count probabilities based on the parts in dna_l.
    
    dna_l   - list of dna instances with mapping_reads defined
    M       - vector of total number of reads at each sample
    
    """
    pseudo_par = np.ones(len(dna_l[0].mapping_reads))
    for dna in dna_l:
        pseudo_par += dna.mapping_reads
    return pseudo_par/(M+1)

def neg_log_probability_l(alpha,dna_l):
    pass

def log_probability(seq,q,M,L):
    """ Calculates the log probability for a time series
    corresponding to one contig and several samples.
    
    seq - dna instance with mapping_reads defined
    q   - vector with relative frequency of genome in each sample
    M   - vector of total number of reads at each sample
    L   - total number of nucleotides in the assembly
    """
    q_prime = q*len(seq)/float(L)
    failed_trials = M-seq.mapping_reads
    return np.sum(seq.mapping_reads*log(q_prime)+
                  failed_trials*log(1-q_prime) + 
                  gammaln(M+1) - 
                  gammaln(seq.mapping_reads+1) - 
                  gammaln(failed_trials+1))
