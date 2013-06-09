#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
import sys
from random import randint
from os import getpid
from itertools import izip
def _clustering(cluster_count, max_iter, run, epsilon, verbose, log_probabilities_func, fit_nonzero_parameters_func, centroids, **kwargs):
    contigs = kwargs["composition"]
    rs = np.random.RandomState(seed=randint(0,10000)+getpid())    
    if not np.any(centroids):
        #centroids = _generate_centroids(contigs,cluster_count, DNA.kmer_hash_count,rs)
        centroids = _generate_kplusplus(contigs, log_probabilities_func,fit_nonzero_parameters_func,cluster_count,DNA.kmer_hash_count,rs)
       
    prev_prob = -np.inf
    prob_diff = np.inf
    iteration = 0

    while (prob_diff >= epsilon and max_iter-iteration > 0):

        cluster_ind = _expectation(contigs, log_probabilities_func, centroids, **kwargs)

        centroids = _maximization(contigs, fit_nonzero_parameters_func, cluster_ind, centroids, rs)
        
        curr_prob = _evaluate_clustering(contigs, log_probabilities_func, cluster_ind, centroids)
        prob_diff = curr_prob - prev_prob 
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        iteration += 1
    #Change back so curr_prob represents the highest probability
    (curr_prob,prev_prob) = (prev_prob,curr_prob)
    if prob_diff < 0:
        print>>sys.stderr, "Kmeans got worse, diff: {0}".format(prob_diff)
    print >> sys.stderr, "Kmeans iterations: {0}".format(iteration)
    clusters = [kwargs["ids"][cluster_ind==i] for i in xrange(len(centroids))]
    return (clusters, curr_prob, centroids)
def _expectation(contigs, log_probabilities_func, centroids, **kwargs):
    prob = log_probabilities_func(contigs,centroids)
    cluster_ind = np.argmax(prob,axis=1)
    return cluster_ind
def _maximization(contigs, fit_nonzero_parameters_func, cluster_ind, centroids,rs):
    centroids[:] = 0
    for i in xrange(len(centroids)):
        #indexes is boolean np.array 
        indexes = cluster_ind==i
        if not indexes.any():
            select_as_centroid = rs.randint(0,len(contigs))
            new_centroid = fit_nonzero_parameters_func([contigs[select_as_centroid]])
            print>>sys.stderr,"cluster {0} was empty in kmeans".format(i)
        else:
            new_centroid = fit_nonzero_parameters_func(contigs[indexes])
        centroids[i] = new_centroid
    return centroids
def _generate_centroids(c_count,c_dim,rs):
    
    centroids = rs.rand(c_count,c_dim)
    centroids /= np.sum(centroids,axis=1,keepdims=True)
    return centroids
def _generate_kplusplus(contigs, log_probabilities_func,fit_nonzero_parameters_func,c_count,c_dim,rs):
    contigs_ind = np.ones(len(contigs)) == 1
    centroids = np.zeros((c_count,c_dim))
    contig_ind = rs.randint(0,len(contigs_ind))
    contigs_ind[contig_ind] = False
    centroids[0] = fit_nonzero_parameters_func(contigs[contig_ind].reshape(1,-1))
    for centroids_index in xrange(1,c_count):
        prob = np.sum(log_probabilities_func(contigs[contigs_ind],centroids[:centroids_index]),axis=1,keepdims=True)
        furthest = np.argmin(prob * np.random.random(prob.shape))
        contigs_ind[furthest] = False
        centroids[centroids_index] = fit_nonzero_parameters_func(contigs[furthest].reshape(1,-1))
    return centroids
def _evaluate_clustering(contigs,log_probabilities_func, cluster_ind, centroids):
    cluster_prob = 0
    for i in xrange(len(centroids)):
        cluster_prob += np.sum(log_probabilities_func(contigs[cluster_ind==i],centroids[i]))
    return cluster_prob


