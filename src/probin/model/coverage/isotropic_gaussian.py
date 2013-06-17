import scipy.stats
import numpy as np

def pdf(mu,contigs=None,**kwargs):
    sigma=kwargs["sigma"]
    qs = np.zeros((len(contigs),len(mu)))
    for i,contig in enumerate(contigs):
        qs[i] = scipy.stats.norm.pdf(contig,loc=mu,scale=sigma).prod(axis=1)
    return qs

#def log_probabilities(mu,contigs=None,**kwargs):
#    return log_pdf(mu,sigma,contigs,**kwargs)

def log_pdf(mu,contigs=None,**kwargs):
    sigma = kwargs["sigma"]
    log_qs = np.zeros((len(contigs),len(mu)))
    for i,contig in enumerate(contigs):
        log_qs[i] = np.log(scipy.stats.norm.pdf(contig,loc=mu,scale=sigma[:len(mu)])).sum(axis=1)
    return log_qs

def fit_parameters(expected_clustering=None,contigs=None,**kwargs):
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
    sigma = kwargs["sigma"]
    if N == 1:
        sigma[:] = 1
    else:
        sigma[:] = 0
        for k in xrange(K):
            sigma[k] = np.array([np.dot(expected_clustering[i,k]*(contigs[i,:]-mu[k,:]),(contigs[i,:]-mu[k,:]).T) for i in xrange(N)]).sum()/n[0,k]

    return mu

