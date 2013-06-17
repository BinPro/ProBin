#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from scipy.special import gammaln
import numpy as np
from random import randint
from os import getpid
import sys


def fit_nonzero_parameters(contigs,expected_clustering=None):
    if expected_clustering == None:
        expected_clustering = np.ones((len(contigs),1))
    pseudo_sig = (expected_clustering.T).dot(contigs) + 1
    pseudo_sig /= np.sum(pseudo_sig,axis=1,keepdims=True)

    return pseudo_sig

def log_probabilities(contigs, p):
    if len(p.shape) == 1:
        p = np.array([p])    
    log_qs = np.zeros((contigs.shape[0],p.shape[0]))
    for i,sign in enumerate(contigs):
        # gammaln produces the natural logarithm of the factorial of i-1
        log_qs[i] = np.sum((sign * np.log(p)) - gammaln(sign+1),axis=1) + gammaln(np.sum(sign)+1)
    return log_qs

def em(contigs, p, K, epsilon, max_iter, **kwargs):
    #N number of contigs
    #D number of features
    #K number of clusters
    N,D = contigs.shape
    
    #initialize with kmeans
    if not np.any(p):
        clustering,_, p = kmeans(contigs, p, K, epsilon, 2, **kwargs)
        
    else:
        print >> sys.stderr, "Not implemented for EM to start with fixed p (centroids)"
        sys.exit(-1)
    
    n = np.array([len(contigs[clustering==c]) for c in xrange(K)],dtype=float)
    
    log_qs = log_probabilities(contigs,p)
    max_log_qs = np.max(log_qs,axis=1,keepdims=True)
    log_qs = np.exp(log_qs - max_log_qs)

    prob_diff = np.inf
    prev_prob = -np.inf
    iteration = 0
    
    while(max_iter - iteration > 0 and prob_diff >= epsilon):
        #================================
        #Expectation
        #================================
        log_qs *= n
        z = log_qs / np.sum(log_qs,axis=1,keepdims=True)
        
        #================================
        #Maximization
        #================================
        p_new = fit_nonzero_parameters(contigs,z)
        
        #================================
        #Evaluation
        #================================
        log_qs = log_probabilities(contigs,p_new)
        max_log_qs = np.max(log_qs,axis=1,keepdims=True)
        log_qs = np.exp(log_qs - max_log_qs)

        curr_prob = np.sum((max_log_qs - np.log(np.sum(z*log_qs,axis=1,keepdims=True))))
        
        n = np.sum(z,axis=0,keepdims=True)

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
    z = log_probabilities(contigs,p)
    #Find each contigs most likely cluster
    clustering = np.argmax(z,axis=1)
    return (clustering, curr_prob, p)

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
    p_new = np.zeros(p.shape)
    prev_prob = -np.inf
    prob_diff = np.inf
    iteration = 0
    
    while (prob_diff >= epsilon and max_iter-iteration > 0):
        #================================
        #Expectation
        #================================
        #Calculate responsibility
        z = log_probabilities(contigs,p)
        #Find each contigs most likely cluster
        clustering_ind = np.argmax(z,axis=1)
        
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
                p_new[K_ind] = fit_nonzero_parameters(contigs[new_centroid])
                print>>sys.stderr,"cluster {0} was empty in kmeans".format(K_ind)
            #Fit the contigs that belong to this cluster to the center
            else:
                p_new[K_ind] = fit_nonzero_parameters(contigs[clustering_K_ind])

        #================================
        #Evaluation
        #================================
        curr_prob = 0
        #for each cluster
        for K_ind in xrange(K):
            #create a boolean array representing that clusters so p[p_ind] is a 2D array (1xD)
            p_ind = np.arange(K) == K_ind
            #calculate log_probabilities of all contigs belonging to cluster k
            curr_prob += np.sum(log_probabilities(contigs[clustering_ind==K_ind],p_new[p_ind]))
        
        prob_diff = curr_prob - prev_prob 
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,p_new) = (p_new,p)
        
        iteration += 1      
    #reverse so curr_prob is the current probability
    (curr_prob,prev_prob) = (prev_prob,curr_prob)
    if prob_diff < 0:
        print>>sys.stderr, "Kmeans got worse, diff: {0}".format(prob_diff)
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        (p,p_new) = (p_new,p)
    print >> sys.stderr, "Kmeans iterations: {0}".format(iteration)
    
    #Get current clustering
    z = log_probabilities(contigs,p)
    #Find each contigs most likely cluster
    clustering = np.argmax(z,axis=1)
    return (clustering, curr_prob, p)
