#!/usr/bin/env python
"""Cluster DNA based on EM algorithm with fixed number of clusters"""
import numpy as np
from itertools import izip
import sys
from multiprocessing import Pool, cpu_count

def cluster(contigs, log_probability_func,fit_nonzero_parameters_func,cluster_count,centroids=None,max_iter=100, repeat=10,epsilon=0.001):
    (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)
    params = [(contigs, log_probability_func,fit_nonzero_parameters_func, cluster_count ,np.copy(centroids), max_iter,epsilon) for _ in xrange(repeat)]
    pool = Pool(processes=cpu_count())
    results = pool.map(_clustering_wrapper, params)
    pool.close()
    return max(results,key=lambda x: x[1])

def _clustering_wrapper(params):
    return _clustering(*params)

def _clustering(contigs, log_probability_func,fit_nonzero_parameters_func, cluster_count ,centroids, max_iter,epsilon):
    if not centroids:
        from probin.binning import kmeans
        (clusters,_,centroids)      = kmeans.cluster(contigs, log_probability_func,fit_nonzero_parameters_func,cluster_count,None,max_iter=3,repeat=1)
        clustering                  = kmeans._expectation(contigs, log_probability_func,centroids)
        expected_cluster_freq       = np.array([len(exp) for exp in clustering],dtype=float)
        clustering_prob, max_log_qs = _evaluate_clustering(centroids, contigs, log_probability_func,expected_cluster_freq)
        
    else:
        print >> sys.stderr, "Not implemented to execute with predefined centroids"
        sys.exit(-1)
    
    it = 0
    while (max_iter-it != 0):

        expected_clustering                = _expectation(contigs,log_probability_func,centroids, expected_cluster_freq, max_log_qs)
        centroids                          = _maximization(contigs, fit_nonzero_parameters_func, centroids, expected_clustering,)
        expected_cluster_freq              = expected_clustering.sum(axis=0,keepdims=True)
        curr_clustering_prob, max_log_qs   = _evaluate_clustering(centroids, contigs, log_probability_func,expected_cluster_freq)
        if (1 - (curr_clustering_prob / clustering_prob) <= epsilon):
            if (curr_clustering_prob < clustering_prob):
                print>>sys.stderr, "EM got worse, previous clustering probability : {0}, current clustering probability: {1}".format( clustering_prob, curr_clustering_prob)
                break
            clustering_prob = curr_clustering_prob
            break
        clustering_prob = curr_clustering_prob
        it += 1
    if not max_iter-it:
        print>>sys.stderr,"EM Finished maximum iteration"
    print >> sys.stderr, "EM executed %s iterations out of %s allowed" % (it,max_iter)
    clusters = [set() for _ in range(cluster_count)]
    [clusters[i].add(contig) for (i,contig) in zip(expected_clustering.argmax(axis=1),contigs)]
    return (clusters, clustering_prob, centroids)

def _expectation(contigs, log_probability_func, centroids,expected_cluster_freq, max_log_qs):
    expected_clustering = np.zeros((len(contigs),len(centroids)))

    for i,contig in enumerate(contigs):
        expected_clustering[i] = [log_probability_func(contig,centroid) for centroid in centroids]
    max_lq = np.max(expected_clustering,axis=1,keepdims=True)
    expected_clustering = np.exp(expected_clustering - max_lq) * expected_cluster_freq
    return expected_clustering / np.sum(expected_clustering,axis=1,keepdims=True)

def _maximization(contigs, fit_nonzero_parameters_func,centroids, expected_clustering):
    for i,exp_cluster in enumerate(expected_clustering.T):
        centroids[i] = fit_nonzero_parameters_func(contigs,expected_clustering=exp_cluster)
    return centroids

def _evaluate_clustering(centroids,contigs, log_probability_func, expected_clustering_freq):
    cluster_prob = 0
    log_qs = np.zero((len(contigs),len(centroids)))
    
    for i,contig in enumerate(contigs):
        log_qs[i] = np.fromiter((log_probability_func(contig,centroid) for centroid in centroids))
    max_log_qs = np.max(log_qs,axis=1,keepdims=True)
        
    exp_diff = np.exp(log_qs - max_log_qs)
    
    
        
#    for (centroid,exp_clust) in izip(centroids,expected_clustering_freq.flatten()):
#        cluster_prob += np.sum(np.array([log_probability_func(contig,centroid) for contig in contigs]) + np.log(exp_clust))
#    return cluster_prob
