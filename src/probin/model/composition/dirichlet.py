"""Implementation of a dirichlet model based on sequence composition"""

from numpy import log
import numpy as np
import sys
from scipy.special import gammaln,psi
from scipy.optimize import fmin_l_bfgs_b, fmin_tnc

def fit_nonzero_parameters(dna_l):
    return np.array(fit_nonzero_parameters_full_output(dna_l)[0])

def fit_nonzero_parameters_full_output(dna_l):
    kmer_hash_count = dna_l[0].kmer_hash_count
    alpha0, pcs = _all_pseudo_counts(dna_l,kmer_hash_count)
    # K is the number of sequences in the cluster
    K,_ = np.shape(pcs)
    # Total number of kmers for each contig
    M = np.sum(pcs,axis=1)
    sys.stderr.write("alpha0: " + str(alpha0) + '\n')
    alpha_bounds = [(0.0,None)]*kmer_hash_count
    alpha_fit = fmin_tnc(neg_log_probability_l,alpha0,fprime=neg_log_probability_l_gradient,args=(pcs,K,M),bounds=alpha_bounds)
    # alpha_fit = fmin_l_bfgs_b(neg_log_probability_l,alpha0,fprime=neg_log_probability_l_gradient,args=(pcs,K,M),bounds=alpha_bounds)
    return alpha_fit

def sophisticated_alpha0(pcs,kmer_hash_count):
    """Estimating the alpha parameters according to a equation given in Minka2012"""
    norms = np.sum(pcs,axis=1)
    pcs_norm = pcs / norms.reshape(-1,1)
    p_expected = np.mean(pcs_norm, axis=0)
    p_var = np.var(pcs_norm, axis=0)
    # The sum of the alphas correspond to the precision of the dirichlet, and is therefore 
    # what makes a difference between the multinomial and the dirichlet-multinomial.
    log_alpha_sum = (1/float(kmer_hash_count-1))*(np.sum(np.log(p_expected*(1-p_expected)/p_var -1))-np.log(p_expected[-1]*(1-p_expected[-1])/p_var[-1] - 1))
    return np.sum(pcs,axis=0)/np.exp(log_alpha_sum)


def neg_log_probability_l_gradient(alpha,pcs,K,M):
    A = np.sum(alpha)
    return -(K*psi(A) - np.sum(psi(A+M)) - K*psi(alpha) + np.sum(psi(pcs+alpha),axis=0))


def neg_log_probability_l(alpha,pcs,K,M):
    """ Negative logarithm of the dirichlet probability 
    
    alpha - vector of dirichlet alpha parameters
    pcs - list of pseudocounts, each corresponding to a sequence
    K - Number of seuquences in the cluster
    M - Total number of kmers for each contig
    """
    A = np.sum(alpha)
    return -(K*gammaln(A) - 
             np.sum(gammaln(A + M)) + 
             np.sum(np.sum(gammaln(pcs+alpha),axis=0) - K*gammaln(alpha)))

def _all_pseudo_counts(dna_l, kmer_hash_count):
    pcs = np.zeros((len(dna_l),kmer_hash_count))
    for index,seq in enumerate(dna_l):
        pcs[index,:] = np.fromiter(seq.pseudo_counts, np.dtype('u4'), count=kmer_hash_count)
    alpha_0 = sophisticated_alpha0(pcs,kmer_hash_count)
    # np.sum(pcs,axis=0)
    return alpha_0,pcs


def log_probability(seq,alpha,pseudo_counts_supplied = False):
    N = np.shape(alpha)[0]
    if pseudo_counts_supplied:
        pc_mat = seq.pseudo_counts_array
    else:
        pc_arr = np.fromiter(seq.pseudo_counts,np.dtype('u4'),count=N)
        pc_mat = pc_arr.reshape((1,N))
    M = np.sum(pc_mat,axis=1)
    return - neg_log_probability_l(alpha,pc_mat,1,M)


def log_probability_test(pseudo_counts,alpha):
    # pseudo_counts given directly mainly for testing purposes
    pcs = np.zeros((1,len(pseudo_counts)))
    pcs[0,:] = np.array(pseudo_counts)
    return - neg_log_probability_l(alpha,pcs)
