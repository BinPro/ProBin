#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
from numpy.random import rand
from collections import Counter

def cluster(contigs, cluster_count,model):
        
    centroids = rand(cluster_count,DNA.kmer_hash_count)
    centroids /= np.sum(centroids,axis=1,keepdims=True)
    clusters = [set() for i in xrange(cluster_count)]
    for contig in contigs:
        classification = _expectation(contig,centroids,model)
        clusters[classification].add(contig)
    
    centroids[::] = _maximization(clusters,model)
    return clusters
    
def _expectation(contig,centroids,model):
    log_prob = {}
    for i,centroid in enumerate(centroids):
        log_prob[model.log_probability(contig.signature,centroid)] = i
    return log_prob[max(log_prob)]

def _maximization(clusters,model):
    centroids_new = np.zeros((len(clusters),DNA.kmer_hash_count))
    for i,cluster in enumerate(clusters):
        clust_sig = Counter()
        [clust_sig.update(contig.signature) for contig in cluster]
        centroids_new[i,:] = model.fit_nonzero_parameters(clust_sig,DNA.kmer_hash_count)
    return centroids_new