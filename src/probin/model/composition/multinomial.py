#!/usr/bin/env python
"""Implementation of a multinomial model based on sequence composition"""
from scipy.special import gammaln
import numpy as np
from collections import Counter
import sys

def em(contigs,p,**kwargs):
    pass

def kmeans(contigs, p, cluster_count, epsilon, max_iter, **kwargs):
    rs = np.random.RandomState(seed=randint(0,10000)+getpid())
    #initialize centroids with random contigs
    if not np.any(p):
        ind = rs.choice(contigs.shape[0],cluster_count,True)
        indx = np.arange(contigs.shape[0])
        p = np.zeros((cluster_count,contigs.shape[1]))
        for i,centroid in enumerate(ind):
            p[i] = fit_nonzero_parameters(contigs[indx==centroid])

    prev_prob = -np.inf
    prob_diff = np.inf
    iteration = 0
    
    while (prob_diff >= epsilon and max_iter-iteration > 0):
        #Expectation
        prob = log_probabilities(contig,centroids)
        clust_ind = np.argmax(prob)
        clusters[clust_ind].add(contig)

        centroids = _maximization(contigs, fit_nonzero_parameters_func, clusters, centroids.shape,rs)
        
        curr_prob = _evaluate_clustering(log_probabilities_func, clusters, centroids)
        prob_diff = curr_prob - prev_prob 
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        iteration += 1

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

    return np.sum((signature_vector * np.log(prob_vector)) - _log_fac(signature_vector)) + _log_fac(np.sum(signature_vector))

def log_probabilities(seq, prob_vectors):
    if len(prob_vectors.shape) == 1:
        prob_vectors = np.array([prob_vectors])
    signature_vector = np.zeros(seq.kmer_hash_count)
    for key,value in seq.signature.iteritems():
        signature_vector[key] = value
    return np.sum((signature_vector * np.log(prob_vectors)) - _log_fac(signature_vector),axis=1) + _log_fac(np.sum(signature_vector))


def _log_fac(i):
    # gammaln produces the natural logarithm of the factorial of i-1
    return gammaln(i+1)
