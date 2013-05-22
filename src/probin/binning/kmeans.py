#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
import sys
from random import randint
from os import getpid
from multiprocessing import Pool, cpu_count


def cluster(contigs, model ,cluster_count,centroids=None,max_iter=100, repeat=10,epsilon=1E-7):    
    (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)
    params = [(contigs, model.log_probabilities, model.fit_nonzero_parameters, cluster_count ,np.copy(centroids), max_iter,epsilon) for _ in xrange(repeat)]
    pool = Pool(processes=cpu_count())
    results = pool.map(_clustering_wrapper, params)
    pool.close()
#    results = [_clustering_wrapper(param) for param in params]
        
    return max(results,key=lambda x: x[1])

def _clustering_wrapper(params):
    return _clustering(*params)

def _clustering(contigs, log_probabilities_func, fit_nonzero_parameters_func, cluster_count ,centroids, max_iter,epsilon):
    rs = np.random.RandomState(seed=randint(0,10000)+getpid())    
    if not np.any(centroids):
       centroids = _generate_kplusplus(contigs, log_probabilities_func,fit_nonzero_parameters_func,cluster_count,DNA.kmer_hash_count,rs)
       
    prev_prob = -np.inf
    prob_diff = np.inf
    iteration = 0

    while (prob_diff > epsilon and max_iter-iteration > 0):

        clusters = _expectation(contigs, log_probabilities_func, centroids)

        centroids = _maximization(contigs, fit_nonzero_parameters_func, clusters, centroids.shape,rs)
        
        curr_prob = _evaluate_clustering(log_probabilities_func, clusters, centroids)
        prob_diff = curr_prob - prev_prob 
        (curr_prob,prev_prob) = (prev_prob,curr_prob)
        iteration += 1
    #Change back so curr_prob represents the highest probability
    (curr_prob,prev_prob) = (prev_prob,curr_prob)
    if prob_diff < 0:
        print>>sys.stderr, "Kmeans got worse, diff: {0}".format(prob_diff)
    print >> sys.stderr, "Kmeans iterations: {0}".format(iteration)
    return (clusters, curr_prob, centroids)

def _expectation(contigs, log_probabilities_func, centroids):
    clusters = [set() for _ in xrange(len(centroids))]
    for contig in contigs:
        prob = log_probabilities_func(contig,centroids)
        clust_ind = np.argmax(prob)
        clusters[clust_ind].add(contig)
    return clusters

def _maximization(contigs, fit_nonzero_parameters_func, clusters, centroids_shape,rs):
    new_centroids = np.zeros(centroids_shape)
    for clust_ind ,clust in enumerate(clusters):
        if not clust:
            select_as_centroid = rs.randint(0,len(contigs))
            new_centroid = fit_nonzero_parameters_func([contigs[select_as_centroid]])
            print>>sys.stderr,"cluster {0} was empty in kmeans".format(clust_ind)
        else:
            new_centroid = fit_nonzero_parameters_func(list(clust))
        new_centroids[clust_ind,:] = new_centroid
    return new_centroids

def _generate_centroids(c_count,c_dim,rs):
    centroids = rs.rand(c_count,c_dim)
    centroids /= np.sum(centroids,axis=1,keepdims=True)
    return centroids

def _generate_kplusplus(contigs, log_probabilities_func,fit_nonzero_parameters_func,c_count,c_dim,rs):
    contigs_ind = range(len(contigs))
    centroids = np.zeros((c_count,c_dim))
    contig_ind = rs.randint(0,len(contigs_ind))
    contigs_ind.remove(contig_ind)
    centroids[0] = fit_nonzero_parameters_func([contigs[contig_ind]])
    for centroids_index in xrange(1,c_count):
        prob = {}
        for contig_ind in contigs_ind:
            sum_prob = sum(log_probabilities_func(contigs[contig_ind],centroids[:centroids_index]) )
            prob[rs.random_sample()*sum_prob] = contig_ind
        furthest = min(prob)
        contig = contigs[prob[furthest]]
        contigs_ind.remove(prob[furthest])
        centroids[centroids_index,:] = fit_nonzero_parameters_func([contig])
    return centroids

def _evaluate_clustering(log_probabilities_func, clusters, centroids):
    cluster_prob = 0
    for i,cluster in enumerate(clusters):
        cluster_prob += sum([log_probabilities_func(contig,centroids[i]) for contig in cluster])
    return cluster_prob


