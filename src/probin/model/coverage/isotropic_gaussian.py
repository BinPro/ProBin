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
    N,L = x.shape
    if expected_clustering is None:
        expected_clustering = np.ones(1,N)
    n = expected_clustering.sum(axis=0)
    K = n.shape[0]
    mu = np.dot(x.T,expected_clustering)
    mu /= n
    sigma = np.zeros(K)
    for k in xrange(K):
        sigma[k] = np.dot(np.array([np.dot(x[i,:]-mu[:,k],(x[i,:]-mu[:,k]).T) for i in xrange(N)]),expected_clustering[:,k])/n[k]

    return mu,sigma
