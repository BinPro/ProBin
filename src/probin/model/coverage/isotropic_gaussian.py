import scipy.stats
import numpy as np

def pdf(x,mu,sigma):
    return scipy.stats.norm.pdf(x,loc=mu,scale=sigma).prod()

def log_probabilities(x,mu_sigma):
    mu,sigma = mu_sigma
    return log_pdf(x,mu,sigma)

def log_pdf(contigs,mu_sigma,**kwargs):
    mu,sigma = mu_sigma
    return np.log(scipy.stats.norm.pdf(contigs,loc=mu,scale=sigma)).sum()

def fit_parameters(contigs,expected_clustering=None,**kwargs):
    # Observe that mu will be 2-dimensional even though 
    # the number of clusters in expected_clustering is 1 or
    # if expected_clustering is undefined.

    N,L = contigs.shape # N: the number of contigs,L: the number of features
    if expected_clustering is None:
         # All contigs are 100% in this single cluster
        expected_clustering = np.ones((N,1))
    n = expected_clustering.sum(axis=0,keepdims=True)
    K = n.shape[1] # The number of clusters
    mu = np.dot(expected_clustering.T,contigs)
    mu /= n.T
    
    if N == 1:
        sigma = np.ones(K)
    else:
        sigma = np.zeros(K)
        for k in xrange(K):
            sigma[k] = np.array([np.dot(expected_clustering[i,k]*(contigs[i,:]-mu[k,:]),(contigs[i,:]-mu[k,:]).T) for i in xrange(N)]).sum()/n[0,k]

    return mu,sigma

