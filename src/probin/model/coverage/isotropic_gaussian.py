import scipy.stats
import numpy as np
import math
import sys

def pdf(x,mu,sigma):
    return scipy.stats.norm.pdf(x,loc=mu,scale=sigma).prod()

def log_probabilities(x,mu,sigma):
    return log_pdf(x,mu,sigma)

def fit_nonzero_parameters(x,**kwargs):
    return fit_parameters(x,**kwargs)

def log_pdf(x,mu,sigma):
    return np.log(scipy.stats.norm.pdf(x,loc=mu,scale=sigma)).sum()

def fit_parameters(x,expected_clustering=None):
    # Observe that mu will be 2-dimensional even though 
    # the number of clusters in expected_clustering is 1 or
    # if expected_clustering is undefined.

    N,L = x.shape # N: the number of contigs,L: the number of features
    if expected_clustering is None:
         # All contigs are 100% in this single cluster
        expected_clustering = np.ones((N,1))
    n = expected_clustering.sum(axis=0,keepdims=True)
    K = n.shape[1] # The number of clusters
    mu = np.dot(expected_clustering.T,x)
    mu /= n.T
    if N == 1:
        sigma = np.ones(K)
    else:
        sigma = np.zeros(K)
        for k in xrange(K):
            sigma[k] = np.array([np.dot(expected_clustering[i,k]*(x[i,:]-mu[k,:]),(x[i,:]-mu[k,:]).T) for i in xrange(N)]).sum()/n[0,k]

    return mu,sigma
