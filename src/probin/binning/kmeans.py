#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
import sys

def cluster(contigs, log_probability_func,fit_nonzero_parameters_func ,cluster_count,centroids=None,max_iter=100, repeat=10,epsilon=0.01):
    (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)    
    for run in xrange(repeat):
        (clusters, clustering_prob, centroids) = _clustering(contigs,  log_probability_func,fit_nonzero_parameters_func, cluster_count ,centroids, max_iter,epsilon)
        (max_clusters, max_clustering_prob,max_centroids) = max([(max_clusters, max_clustering_prob, max_centroids), (clusters, clustering_prob, centroids)],key=lambda x: x[1])
    return (max_clusters, max_clustering_prob, max_centroids)
def _clustering(contigs,  log_probability_func,fit_nonzero_parameters_func, cluster_count ,centroids, max_iter,epsilon):
    if centroids is None:
       centroids = _generate_kplusplus(contigs, log_probability_func,fit_nonzero_parameters_func,cluster_count,DNA.kmer_hash_count)
    clustering_prob = -np.inf
    cluster_different = True
    
    while (cluster_different and max_iter != 0):

        clusters = _expectation(contigs, log_probability_func, centroids)

        centroids = _maximization(contigs, fit_nonzero_parameters_func, clusters, centroids.shape)
        
        curr_clustering_prob = _evaluate_clustering(log_probability_func, clusters, centroids)
        
        if (1-curr_clustering_prob/clustering_prob <= epsilon):
            cluster_different = False
            if (curr_clustering_prob < clustering_prob):
                print>>sys.stderr, "Kmeans got worse, previous clustering probability : {0}, current clustering probability: {1}".format( clustering_prob, curr_clustering_prob)
        clustering_prob = curr_clustering_prob
        max_iter -= 1
    if not max_iter:
        print>>sys.stderr,"Kmeans Finished maximum iteration"
    return (clusters, clustering_prob, centroids)

def _expectation(contigs, log_probability_func, centroids):
    clusters = [set() for _ in xrange(len(centroids))]
    for contig in contigs:
        prob = [log_probability_func(contig,centroid) for centroid in centroids]
        clust_ind = np.argmax(prob)
        clusters[clust_ind].add(contig)
    return clusters

def _maximization(contigs, fit_nonzero_parameters_func, clusters, centroids_shape):
    new_centroids = np.zeros(centroids_shape)
    for clust_ind ,clust in enumerate(clusters):
        if not clust:
            select_as_centroid = np.random.randint(0,len(contigs))
            new_centroid = fit_nonzero_parameters_func([contigs[select_as_centroid]])
            print>>sys.stderr,"cluster {0} was empty in kmeans".format(clust_ind)
        else:
            new_centroid = fit_nonzero_parameters_func(list(clust))
        new_centroids[clust_ind,:] = new_centroid
    return new_centroids

def _generate_centroids(c_count,c_dim):
    centroids = np.random.rand(c_count,c_dim)
    centroids /= np.sum(centroids,axis=1,keepdims=True)
    return centroids

def _generate_kplusplus(contigs, log_probability_func,fit_nonzero_parameters_func,c_count,c_dim):
    contigs_ind = range(len(contigs))
    centroids = np.zeros((c_count,c_dim))
    contig_ind = np.random.randint(0,len(contigs_ind))
    contigs_ind.remove(contig_ind)
    centroids[0,:] = fit_nonzero_parameters_func([contigs[contig_ind]])
    for centroids_index in xrange(1,c_count):
        prob = {}
        for contig_ind in contigs_ind:
            sum_prob = sum([log_probability_func(contigs[contig_ind],centroid) for centroid in centroids[:centroids_index]])
            prob[np.random.random()*sum_prob] = contig_ind
        furthest = min(prob)
        contig = contigs[prob[furthest]]
        contigs_ind.remove(prob[furthest])
        centroids[centroids_index,:] = fit_nonzero_parameters_func([contig])
    return centroids

def _evaluate_clustering(log_probability_func, clusters, centroids):
    cluster_prob = 0
    for i,cluster in enumerate(clusters):
        cluster_prob += sum([log_probability_func(contig,centroids[i]) for contig in cluster])
    return cluster_prob


