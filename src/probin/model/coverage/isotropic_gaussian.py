import scipy.stats
import numpy as np
from random import randint
from os import getpid
import sys

def pdf(contigs,p,sigma):
    qs = np.zeros((contigs.shape[0],p.shape[0]))
    for i,contig in enumerate(contigs):
        qs[i] = scipy.stats.norm.pdf(contig,loc=p,scale=sigma).prod(axis=1)
    return qs

#def log_probabilities(mu,contigs=None,**kwargs):
#    return log_pdf(mu,sigma,contigs,**kwargs)

#def fit_nonzero_parameters(x,**kwargs):
#    return fit_parameters(x,**kwargs)

def log_pdf(contigs,p,sigma):
    log_qs = np.zeros((contigs.shape[0],p.shape[0]))
    for i,contig in enumerate(contigs):
        log_qs[i] = np.log(scipy.stats.norm.pdf(contig,loc=p,scale=sigma)).sum(axis=1)
    return log_qs

def fit_parameters(contigs,expected_clustering=None):
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
        sigma = np.ones((K,1))
    else:
        #sigma = np.zeros((K,1))
        sigma = np.ones((K,1))
        #for k in xrange(K):
        #    sigma[k] = np.array([np.dot(expected_clustering[i,k]*(contigs[i,:]-mu[k,:]),(contigs[i,:]-mu[k,:]).T) for i in xrange(N)]).sum()/n[0,k]

    return mu,sigma


def em(contigs, p, K, epsilon, max_iter, **kwargs):
    #N number of contigs
    #D number of features
    #K number of clusters
    N,D = contigs.shape
    
    #initialize with kmeans
    if not np.any(p):
        clustering,_, (p,sigma) = kmeans(contigs, p, K, epsilon, 2, **kwargs)
        
    else:
        print >> sys.stderr, "Not implemented for EM to start with fixed p (centroids)"
        sys.exit(-1)
    
    n = np.array([len(contigs[clustering==c]) for c in xrange(K)],dtype=float)

    #Calculate z ones. it is then calculated in the evaluation and reused in next iteration
    log_qs = log_pdf(contigs,p,sigma)
    max_log_qs = np.max(log_qs,axis=1,keepdims=True)
    log_qs = np.exp(log_qs - max_log_qs)

    #initialize
    prob_diff = np.inf
    prev_prob = -np.inf
    iteration = 0
    
    while(max_iter - iteration > 0 and prob_diff >= epsilon):
        #================================
        #Expectation
        #================================
        #Nothing since z values are calculated in evaluation in previous iteration
        z = (log_qs*n) / np.sum((log_qs*n),axis=1,keepdims=True)
        
        #================================
        #Maximization
        #================================
        p_new, sigma = fit_parameters(contigs,z)
        n = np.sum(z,axis=0,keepdims=True)

        #================================
        #Evaluation
        #================================
        log_qs = log_pdf(contigs,p_new,sigma)
        max_log_qs = np.max(log_qs,axis=1,keepdims=True)
        log_qs = np.exp(log_qs - max_log_qs)

        curr_prob = np.sum((max_log_qs - np.log(np.sum(z*log_qs,axis=1,keepdims=True))))
        
        prob_diff = curr_prob - prev_prob

        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,p_new) = (p_new,p)
        iteration += 1

    (curr_prob,prev_prob) = (prev_prob,curr_prob)
    print >> sys.stderr, "EM iterations: {0}, difference: {1}".format(iteration, prob_diff)
    if prob_diff < 0:
        print >> sys.stderr, "EM got worse, diff: {0}".format(prob_diff)
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,p_new) = (p_new,p)
        
    #Get current clustering
    log_qs = log_pdf(contigs,p,sigma)
    #Find each contigs most likely cluster
    clustering = np.argmax(log_qs,axis=1)
    return (clustering, curr_prob, (p,sigma))

def kmeans(contigs, p, K, epsilon, max_iter, **kwargs):
    rs = np.random.RandomState(seed=randint(0,10000)+getpid())
    #N number of contigs
    #D number of features
    #K number of clusters
    N,D = contigs.shape
    #initialize centroids with random contigs
    sigma = np.zeros((K,1))
    if not np.any(p):
        ind = rs.choice(N,K,True)
        indx = np.arange(N)
        p = np.zeros((K,D))
        for i,centroid in enumerate(ind):
            p[i],sigma[i] = fit_parameters(contigs[indx==centroid])
            
    #Calculate responsibility, I don't normalize since 
    #I am only picking the greatest one and using that.
    log_qs = log_pdf(contigs,p,sigma)
    #Find each contigs most likely cluster
    clustering_ind = np.argmax(log_qs,axis=1)
    
    #Initialize
    p_new = np.zeros(p.shape)
    prev_prob = -np.inf
    prob_diff = np.inf
    iteration = 0
    
    while (prob_diff >= epsilon and max_iter-iteration > 0):
        #================================
        #Expectation
        #================================
        #Nothing done since the log_qs and clustering_ind are
        #already calculated in previous iteration
        
        #================================
        #Maximization
        #================================
        # For ecah cluster
        for K_ind in xrange(K):
            #Gives boolean array with true for contigs belonging to cluster K
            clustering_K_ind = clustering_ind == K_ind
            
            #Fit the contigs that belong to this cluster to the center
            if clustering_K_ind.any():
                p_new[K_ind],sigma[K_ind] = fit_parameters(contigs[clustering_K_ind])
            #if empty, pick random contig to represent that clusters
            else:
                new_centroid = np.arange(N) == rs.randint(0,N)
                p_new[K_ind],sigma[K_ind] = fit_parameters(contigs[new_centroid])
                print>>sys.stderr,"cluster {0} was empty in kmeans".format(K_ind)

        #================================
        #Evaluation
        #================================
        curr_prob = 0
        
        #Calculate responsibility, I don't normalize since 
        #I am only picking the greatest one and using that.
        #These are the same number we would get from the expectation step in
        #next iteration so they are reused there.
        log_qs = log_pdf(contigs,p_new,sigma)
        clustering_ind = np.argmax(log_qs,axis=1)
        
        #for each cluster
        for K_ind in xrange(K):
            #create a boolean array representing that clusters so p[p_ind] is a 2D array (1xD)
            p_ind = np.arange(K) == K_ind
            #calculate log_probabilities of all contigs belonging to cluster k
            curr_prob += np.sum(log_pdf(contigs[clustering_ind==K_ind],p_new[p_ind],sigma[p_ind]))
        
        prob_diff = curr_prob - prev_prob 
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,p_new) = (p_new,p)        
        iteration += 1    
        print iteration
        
        
    #reverse so curr_prob is the current probability
    (curr_prob,prev_prob) = (prev_prob,curr_prob)
    if prob_diff < 0:
        print>>sys.stderr, "Kmeans got worse, diff: {0}".format(prob_diff)
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,p_new) = (p_new,p)
    print >> sys.stderr, "Kmeans iterations: {0}".format(iteration)
    
    #Get current clustering, not normalized since we pick the max one
    log_qs = log_pdf(contigs,p,sigma)
    #Find each contigs most likely cluster
    clustering = np.argmax(log_qs,axis=1)
    
    return (clustering, curr_prob, (p, sigma))