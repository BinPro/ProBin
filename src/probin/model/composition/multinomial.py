#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from scipy.special import gammaln
import numpy as np
from collections import Counter
from random import randint
from os import getpid
import sys

def em(contigs,p,**kwargs):
    pass

def kmeans(contigs, p, K, epsilon, max_iter, **kwargs):
    rs = np.random.RandomState(seed=randint(0,10000)+getpid())
    #N number of contigs
    #D number of features
    #K number of clusters
    N,D = contigs.shape
    #initialize centroids with random contigs
    if not np.any(p):
        ind = rs.choice(N,K,True)
        indx = np.arange(N)
        p = np.zeros((K,D))
        for i,centroid in enumerate(ind):
            p[i] = fit_nonzero_parameters(contigs[indx==centroid])
            
    new_p = np.zeros(p.shape)
    z = np.zeros((N,K))
    prev_prob = -np.inf
    prob_diff = np.inf
    iteration = 0
    
    while (prob_diff >= epsilon and max_iter-iteration > 0):
        #================================
        #Expectation
        #================================
        #Calculate responsibility
        z[:] = log_probabilities(contigs,p)
        #Find each contigs most likely cluster
        clustering_ind = np.argmax(z)
        
        #================================
        #Maximization
        #================================
        # For ecah cluster
        for K_ind in xrange(K):
            #Gives boolean array with true for contigs belonging to cluster K
            clustering_K_ind = clustering_ind == K_ind
            #if empty, pick random contig to represent that clusters
            if not clustering_K_ind.any():
                new_centroid = np.arange(N) == rs.randint(0,N)
                new_p[K_ind] = fit_nonzero_parameters(contigs[new_centroid])
                print>>sys.stderr,"cluster {0} was empty in kmeans".format(K_ind)
            #Fit the contigs that belong to this cluster to the center
            else:
                new_p[K_ind] = fit_nonzero_parameters(contigs[clustering_K_ind])

        #================================
        #Evaluation
        #================================
        curr_prob = 0
        #for each cluster
        for K_ind in xrange(K):
            #create a boolean array representing that clusters so p[p_ind] is a 2D array (1xD)
            p_ind = np.arange(K) == K_ind
            #calculate log_probabilities of all contigs belonging to cluster k
            curr_prob += np.sum(log_probabilities(contigs[clustering_ind==K_ind],new_p[p_ind]))
        
        prob_diff = curr_prob - prev_prob 
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,new_p) = (new_p,p)
        
        iteration += 1        
    
    if prob_diff < 0:
        print>>sys.stderr, "Kmeans got worse, diff: {0}".format(prob_diff)
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,new_p) = (new_p,p)
    print >> sys.stderr, "Kmeans iterations: {0}".format(iteration)
    
    #Get current clustering
    z[:] = log_probabilities(contigs,p)
    #Find each contigs most likely cluster
    clustering = np.argmax(z)
    return (clustering, curr_prob, p)

def fit_parameters(dna_l):
    sig = Counter()
    [sig.update(part.signature) for part in dna_l]
    par = np.zeros(dna_l[0].kmer_hash_count)
    for key,cnt in sig.iteritems():
        par[key] += cnt
    par /= np.sum(par)
    return par

def fit_nonzero_parameters(dna_l,expected_clustering=None):
    pseudo_sig = np.zeros((len(dna_l),dna_l[0].kmer_hash_count))
    for i,dna in enumerate(dna_l):
        pseudo_sig[i,:] = np.fromiter(dna.pseudo_counts,dtype=np.int) - 1
    if expected_clustering == None:
        expected_clustering = np.ones((1,len(dna_l)))
    else:
        expected_clustering= expected_clustering.T
    pseudo_sig = expected_clustering.dot(pseudo_sig)
    pseudo_sig += 1
    pseudo_sig /= np.sum(pseudo_sig,axis=1,keepdims=True)
    if len(pseudo_sig) == 1:
        pseudo_sig = pseudo_sig[0]
    return pseudo_sig

def log_probability(seq, prob_vector):
    signature_vector = np.zeros(np.shape(prob_vector))

    for key,value in seq.signature.iteritems():
        signature_vector[key] = value

    return np.sum((signature_vector * np.log(prob_vector)) - gammaln(signature_vector+1)) + gammaln(np.sum(signature_vector)+1)

def log_probabilities(contigs, prob_vectors):
    if len(prob_vectors.shape) == 1:
        prob_vectors = np.array([prob_vectors])    
    log_qs = np.zeros((contigs.shape[0],prob_vectors.shape[0]))
    for i,sign in enumerate(contigs):
        # gammaln produces the natural logarithm of the factorial of i-1
        log_qs[i] = np.sum((sign * np.log(prob_vectors)) - gammaln(sign+1),axis=1) + gammaln(np.sum(sign)+1)
    return log_qs

def _log_fac(i):
    # gammaln produces the natural logarithm of the factorial of i-1
    return gammaln(i+1)
