#!/usr/bin/env python
"""Cluster DNA based on EM algorithm with fixed number of clusters"""
import numpy as np
from itertools import izip
import sys

def cluster(contigs,model,cluster_count,centroids=None,max_iter=100, repeat=10,epsilon=0.01):
    (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)
    
    for run in xrange(repeat):
        (clusters, clustering_prob, new_centroids) = _clustering(contigs, model, cluster_count ,centroids, max_iter,epsilon)
        (max_clusters, max_clustering_prob,max_centroids) = max([(max_clusters, max_clustering_prob, max_centroids), (clusters, clustering_prob, new_centroids)],key=lambda x: x[1])
    return (max_clusters, max_clustering_prob, max_centroids)

def _clustering(contigs, model, cluster_count ,centroids, max_iter,epsilon):
    if not centroids:
        from probin.binning import kmeans
        (clusters,_,centroids) = kmeans.cluster(contigs,model,cluster_count,None,max_iter=3,repeat=2)
        expected_cluster_freq = np.ones((1,cluster_count))/cluster_count
    else:
        clusters = _expectation(contigs,centroids)
    clustering_prob = -np.inf
    while (max_iter != 0):

        expected_clustering     = _expectation(contigs,model,centroids, expected_cluster_freq)
        centroids               = _maximization(contigs, model, centroids, expected_clustering,)
        expected_cluster_freq   = expected_clustering.sum(axis=0,keepdims=True)
        curr_clustering_prob    = _evaluate_clustering(centroids, contigs, model,expected_cluster_freq)
        
        if (1 - curr_clustering_prob / clustering_prob <= epsilon):
            if (curr_clustering_prob < clustering_prob):
                print>>sys.stderr, "EM got worse, previous clustering probability : {0}, current clustering probability: {1}".format( clustering_prob, curr_clustering_prob)
                break
            clustering_prob = curr_clustering_prob
            break
        clustering_prob = curr_clustering_prob
        max_iter -= 1
    if not max_iter:
        print>>sys.stderr,"EM Finished maximum iteration"
    clusters = [set() for _ in range(cluster_count)]
    [clusters[i].add(contig) for (i,contig) in zip(expected_clustering.argmax(axis=1),contigs)]
    return (clusters, clustering_prob, centroids)

def _expectation(contigs, model, centroids,expected_cluster_freq):
    expected_clustering = np.zeros((len(contigs),len(centroids)))

    for i,contig in enumerate(contigs):
        expected_clustering[i] = [model.log_probability(contig,centroid) for centroid in centroids]
    max_lq = np.max(expected_clustering,axis=1,keepdims=True)
    expected_clustering = np.exp(expected_clustering - max_lq) * expected_cluster_freq
    return expected_clustering / np.sum(expected_clustering,axis=1,keepdims=True)

def _maximization(contigs, model,centroids, expected_clustering):
    for i,exp_cluster in enumerate(expected_clustering.T):
        centroids[i] = model.fit_nonzero_parameters(contigs,expected_clustering=exp_cluster)
    return centroids

def _evaluate_clustering(centroids,contigs, model, expected_clustering_freq):
    cluster_prob = 0
    for (centroid,exp_clust) in izip(centroids,expected_clustering_freq.T):
        cluster_prob += np.sum(np.array([model.log_probability(contig,centroid) for contig in contigs]) + np.log(exp_clust))
    return cluster_prob
