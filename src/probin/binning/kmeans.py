#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
from collections import Counter

def cluster(contigs, cluster_count,model,centroids=None):
    if centroids is None:
       centroids = _generate_centroids(cluster_count,DNA.kmer_hash_count)
    cluster_different = True
    while (cluster_different):
        clusters = [set() for _ in xrange(cluster_count)]

        #Expectations
        for contig in contigs:
            prob = [model.log_probability(contig.signature,centroid) for centroid in centroids]
            clust_ind = np.argmax(prob)
            clusters[clust_ind].add(contig)
        
        #Maximization
        new_centroids = np.zeros(centroids.shape)
        for clust_ind ,clust in enumerate(clusters):
            new_centroid = Counter()
            [new_centroid.update(contig.signature) for contig in clust]
            new_centroids[clust_ind,:] = model.fit_nonzero_parameters(new_centroid,DNA.kmer_hash_count)
        if (new_centroids == centroids).all():
            cluster_different = False
        (new_centroids,centroids) = (centroids,new_centroids)
    return centroids
    
def _generate_centroids(c_count,c_dim):
    centroids = np.random.rand(c_count,c_dim)
    centroids /= np.sum(centroids,axis=1,keepdims=True)
    return centroids




