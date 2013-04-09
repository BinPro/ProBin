"""Implementation of a dirichlet model based on sequence composition"""

from numpy import log
import numpy as np
from scipy.special import gammaln
from scipy.optimize import fmin_l_bfgs_b

def fit_nonzero_parameters(dna_l):
    return np.array(fit_nonzero_parameters_full_output(dna_l)[0])

def fit_nonzero_parameters_full_output(dna_l):
    kmer_hash_count = dna_l[0].kmer_hash_count
    alpha0, pcs = _all_pseudo_counts(dna_l,kmer_hash_count)
    alpha_bounds = [(0.0,None)]*kmer_hash_count
    alpha_fit = fmin_l_bfgs_b(neg_log_probability_l,alpha0,args=(pcs,),bounds=alpha_bounds, approx_grad=True, epsilon=1e-12)
    return alpha_fit

def neg_log_probability_l(alpha,pcs):
    A = np.sum(alpha)

    # N is the number of sequences in the sample
    N,_ = np.shape(pcs) # Different meaning than before

    # Total number of kmers for each contig
    M = np.sum(pcs,axis=1)

    return -(N*gammaln(A) - 
             np.sum(gammaln(A + M)) + 
             np.sum(np.sum(gammaln(pcs+alpha),axis=0) - N*gammaln(alpha)))

def _all_pseudo_counts(dna_l, kmer_hash_count):
    pcs = np.zeros((len(dna_l),kmer_hash_count))
    for index,seq in enumerate(dna_l):
        pcs[index,:] = np.fromiter(seq.pseudo_counts, np.dtype('u4'), count=kmer_hash_count)
    alpha_0 = np.sum(pcs,axis=0)
    return alpha_0,pcs

def log_probability(seq,alpha):
    N = np.shape(alpha)[0]
    pc_arr = np.fromiter(seq.pseudo_counts,np.dtype('u4'),count=N)
    pc_mat = pc_arr.reshape((1,N))
    return - neg_log_probability_l(alpha,pc_mat)

def log_probability_test(pseudo_counts,alpha):
    # pseudo_counts given directly mainly for testing purposes
    pcs = np.zeros((1,len(pseudo_counts)))
    pcs[0,:] = np.array(pseudo_counts)
    return - neg_log_probability_l(alpha,pcs)
