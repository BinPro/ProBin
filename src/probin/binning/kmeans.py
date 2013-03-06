#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
from numpy.random import rand
def cluster(contigs, cluster_count,model):
        
    centroids = rand(cluster_count,DNA.kmer_hash_count)
    centroids /= np.sum(centroids,axis=1,keepdims=True)
    clusters = [set() for i in xrange(cluster_count)]
    for contig in contigs:
        classification = _expectation(contig,centroids,model)
        clusters[classification].add(contig)
    
    _maximization(centroids,clusters)
    
    return clusters
    
def _expectation(contig,centroids,model):
    log_prob = {}
    for i,centroid in enumerate(centroids):
        log_prob[model.log_probability(contig.signature,centroid)] = i
    return log_prob[max(log_prob)]

def _maximization(centroids,clusters):
    
    