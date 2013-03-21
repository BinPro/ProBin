#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
from collections import Counter
import sys

def cluster(contigs, model, cluster_count ,centroids=None, max_iter=100, repeat=10):
    (max_clustering_prob,max_centroids,max_clusters) = (-np.inf,None,None)
    if repeat != 1:
        (max_clustering_prob,max_centroids,max_clusters) = cluster(contigs, model, cluster_count ,centroids=None, max_iter=100, repeat=repeat-1)
    if centroids is None:
       centroids = _generate_kplusplus(contigs,model,cluster_count,DNA.kmer_hash_count)
    clustering_prob = -np.inf
    cluster_different = True
    
    while (cluster_different and max_iter != 0):

        clusters = _expectation(contigs,model,centroids)

        centroids = _maximization(contigs, model, clusters, centroids.shape)
        
        curr_clustering_prob = _evaluate_clustering(centroids, clusters, model)
        
        if (curr_clustering_prob <= clustering_prob):
            cluster_different = False
            if (curr_clustering_prob < clustering_prob):
                print>>sys.stderr, "Clustering got worse, previous clustering probability : {0}, current clustering probability: {1}".format( clustering_prob, curr_clustering_prob)
        clustering_prob = curr_clustering_prob
        max_iter -= 1
    (curr_max_clust_prob, curr_max_centr, curr_max_clust) = max([(max_clustering_prob, max_centroids, max_clusters), (clustering_prob, centroids, clusters)],key=lambda x: x[0])
    return (curr_max_clust_prob, curr_max_centr, curr_max_clust)

def _expectation(contigs, model, centroids):
    clusters = [set() for _ in xrange(len(centroids))]
    for contig in contigs:
        prob = [model.log_probability(contig.signature,centroid) for centroid in centroids]
        clust_ind = np.argmax(prob)
        clusters[clust_ind].add(contig)
    return clusters
    
def _maximization(contigs, model, clusters, centroids_shape):
    new_centroids = np.zeros(centroids_shape)
    for clust_ind ,clust in enumerate(clusters):
        if not clust:
            select_as_centroid = np.random.randint(0,len(contigs))
            new_centroid = model.fit_nonzero_parameters(contigs[select_as_centroid].signature,DNA.kmer_hash_count)
        else:
            new_centroid_count = Counter()
            [new_centroid_count.update(contig.signature) for contig in clust]
            new_centroid = model.fit_nonzero_parameters(new_centroid_count,DNA.kmer_hash_count)
        new_centroids[clust_ind,:] = new_centroid
    return new_centroids

def _generate_centroids(c_count,c_dim):
    centroids = np.random.rand(c_count,c_dim)
    centroids /= np.sum(centroids,axis=1,keepdims=True)
    return centroids

def _generate_kplusplus(contigs,model,c_count,c_dim):
    contigs_ind = range(len(contigs))
    centroids = np.zeros((c_count,c_dim))
    contig_ind = np.random.randint(0,len(contigs_ind))
    contigs_ind.remove(contig_ind)
    centroids[0,:] = model.fit_nonzero_parameters(contigs[contig_ind].signature,DNA.kmer_hash_count)
    for centroids_index in xrange(1,c_count):
        prob = {}
        for contig_ind in contigs_ind:
            sum_prob = sum([model.log_probability(contigs[contig_ind].signature,centroid) for centroid in centroids[:centroids_index]])
            prob[np.random.random()*sum_prob] = contig_ind
        furthest = min(prob)
        contig = contigs[prob[furthest]]
        contigs_ind.remove(prob[furthest])
        centroids[centroids_index,:] = model.fit_nonzero_parameters(contig.signature,DNA.kmer_hash_count)
    return centroids
    
def _evaluate_clustering(centroids,clusters, model):
    cluster_prob = 0
    for i,cluster in enumerate(clusters):
        cluster_prob += sum([model.log_probability(contig.signature,centroids[i]) for contig in cluster])
    return cluster_prob


