#!/usr/bin/env python
"""Cluster DNA based on EM algorithm with fixed number of clusters"""
import numpy as np
from itertools import izip
import sys

def cluster(contigs,model,cluster_count,centroids=None,max_iter=100, repeat=10):
    (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)
    
    for run in xrange(repeat):
        (clusters, clustering_prob, new_centroids) = _clustering(contigs, model, cluster_count ,centroids, max_iter,repeat)
        (max_clusters, max_clustering_prob,max_centroids) = max([(max_clusters, max_clustering_prob, max_centroids), (clusters, clustering_prob, new_centroids)],key=lambda x: x[1])
    return (max_clusters, max_clustering_prob, max_centroids)
    
def _clustering(contigs, model, cluster_count ,centroids, max_iter,repeat):
    if not centroids:
        from probin.binning import kmeans
        (clusters,clustering_prob,centroids) = kmeans.cluster(contigs,model,cluster_count,centroids,max_iter,repeat)
    else:
        clusters = _expectation(contigs,centroids)
        
    cluster_different = True
    
    while (cluster_different and max_iter != 0):
        
        cluster_freq            = [len(cluster) for cluster in clusters]
        clusters                = _expectation(contigs,model,centroids, cluster_freq)
        centroids               = _maximization(contigs, model, clusters, centroids.shape)
        curr_clustering_prob    = _evaluate_clustering(centroids, clusters, model)
        
        if (curr_clustering_prob <= clustering_prob):
            cluster_different = False
            if (curr_clustering_prob < clustering_prob):
                print>>sys.stderr, "Clustering got worse, previous clustering probability : {0}, current clustering probability: {1}".format( clustering_prob, curr_clustering_prob)
        clustering_prob = curr_clustering_prob
        max_iter -= 1
    if not max_iter:
        print>>sys.stderr,"Finished maximum iteration"
    return (clusters, clustering_prob, centroids)

def _expectation(contigs, model, centroids,cluster_freq):
    clusters = [set() for _ in xrange(len(centroids))]
    for contig in contigs:
        # Since we use log probability we need to divide by the cluster frequency rather than multiply by it.
        prob = [model.log_probability(contig,centroid)/clust_freq for (centroid,clust_freq) in izip(centroids,cluster_freq)]
        clust_ind = np.argmax(prob)
        clusters[clust_ind].add(contig)
    return clusters

def _maximization(contigs, model, clusters, centroids_shape):
    new_centroids = np.zeros(centroids_shape)
    for clust_ind ,clust in enumerate(clusters):
        if not clust:
            select_as_centroid = np.random.randint(0,len(contigs))
            new_centroid = model.fit_nonzero_parameters([contigs[select_as_centroid]])
        else:
            new_centroid = model.fit_nonzero_parameters(list(clust))
        new_centroids[clust_ind,:] = new_centroid
    return new_centroids


def _evaluate_clustering(centroids,clusters, model):
    cluster_prob = 0
    for i,cluster in enumerate(clusters):
        cluster_prob += sum([model.log_probability(contig,centroids[i]) for contig in cluster])
    return cluster_prob


