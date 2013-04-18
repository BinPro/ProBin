"""Implementation of a dirichlet model based on sequence composition"""

from numpy import log
import numpy as np
import sys
from scipy.special import gammaln,psi
from scipy.optimize import fmin_l_bfgs_b, fmin_tnc

def fit_nonzero_parameters(dna_l,algorithm="tnc"):
    return np.array(fit_nonzero_parameters_full_output(dna_l,algorithm=algorithm)[0])

def fit_nonzero_parameters_full_output(dna_l, algorithm="tnc"):
    """ Approximate maximum likelihood estimates for dirichlet.

    dna_l   -- list of probin.dna.DNA instances

    Keyword Arguments:
    algorithm -- choose numerical optimization algorithm, e.g. "bfgs"

    """
    kmer_hash_count = dna_l[0].kmer_hash_count
    alpha0, pcs = _all_pseudo_counts(dna_l,kmer_hash_count)
    # K is the number of sequences in the cluster
    K,_ = np.shape(pcs)
    # Total number of kmers for each contig
    M = np.sum(pcs,axis=1)
    alpha_bounds = [(0.0,None)]*kmer_hash_count
    if algorithm == "bfgs":
        alpha_fit = fmin_l_bfgs_b(neg_log_probability_l,alpha0,fprime=neg_log_probability_l_gradient,args=(pcs,K,M),bounds=alpha_bounds)
    else:
        alpha_fit = fmin_tnc(neg_log_probability_l,alpha0,fprime=neg_log_probability_l_gradient,args=(pcs,K,M),bounds=alpha_bounds)    
    return alpha_fit

def sophisticated_alpha0(pcs,kmer_hash_count):
    """Estimating alpha parameters according to an equation given in Minka2012.

    pcs             -- pseudo counts for each kmer and each sequence,
                       at least two sequence is needed.
    kmer_hash_count -- number of unique kmers

    """
    norms = np.sum(pcs,axis=1)
    pcs_norm = pcs / norms.reshape(-1,1)
    p_expected = np.mean(pcs_norm, axis=0)
    p_var = np.var(pcs_norm, axis=0)
    if np.min(p_var) == 0.0:
        p_var += np.max(p_var)**2
    # The sum of the alphas correspond to the precision of the dirichlet, and is therefore 
    # what makes a difference between the multinomial and the dirichlet-multinomial.
    log_alpha_sum = (1/float(kmer_hash_count-1))*(np.sum(np.log(p_expected*(1-p_expected)/p_var -1))-np.log(p_expected[-1]*(1-p_expected[-1])/p_var[-1] - 1))
    return np.sum(pcs,axis=0)/np.exp(log_alpha_sum)


def neg_log_probability_l_gradient(alpha,sigs,K,M):
    """ Gradient for the negative logarithm of the dirichlet probability.

    alpha -- vector of dirichlet alpha parameters
    sigs  -- genomic signatures, each corresponding to a sequence
    K     -- Number of seuquences in the cluster
    M     -- Total number of kmers for each contig

    """
    A = np.sum(alpha)
    return -(K*psi(A) - np.sum(psi(A+M)) - K*psi(alpha) + np.sum(psi(sigs+alpha),axis=0))


def neg_log_probability_l(alpha,sigs,K,M):
    """ 
    Negative logarithm of the dirichlet probability 
    
    alpha -- vector of dirichlet alpha parameters
    sigs  -- genomic signatures, each corresponding to a sequence
    K     -- Number of seuquences in the cluster
    M     -- Total number of kmers for each contig

    """

    A = np.sum(alpha)
    return -(K*gammaln(A) - 
             np.sum(gammaln(A + M)) + 
             np.sum(np.sum(gammaln(sigs+alpha),axis=0) - K*gammaln(alpha)))

def _all_pseudo_counts(dna_l, kmer_hash_count):
    """
    Calculates initial guess for alpha and pseudo counts for dna_l

    dna_l - list of dna sequences
    kmer_hash_count - number of unique kmers

    """
    pcs = np.zeros((len(dna_l),kmer_hash_count))
    for index,seq in enumerate(dna_l):
        pcs[index,:] = np.fromiter(seq.pseudo_counts, np.dtype('u4'), count=kmer_hash_count)
    alpha_0 = sophisticated_alpha0(pcs,kmer_hash_count)
    # np.sum(pcs,axis=0)
    return alpha_0,pcs


def log_probability(seq,alpha):
    """ Calculates the logarithm of the dirichlet probability.

    seq    -- probin.dna.DNA instance.
    alpha  -- numpy vector with dirichlet parameters.

    """
    N = np.shape(alpha)[0]

    sig_arr = np.zeros((1,N))
    for key,cnt in seq.signature.iteritems():
        sig_arr[0,key] += cnt

    M = np.sum(sig_arr,axis=1)
    return - neg_log_probability_l(alpha,sig_arr,1,M)


def log_probability_test(signature,alpha):
    """ Calculate probability of signature directly. 

    Mainly for testing purposes

    """
    sig = signature.reshape((1,len(signature)))
    return -neg_log_probability_l(alpha,sig,1,np.sum(signature))
