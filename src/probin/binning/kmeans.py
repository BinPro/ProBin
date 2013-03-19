#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
from collections import Counter
import sys

def cluster(contigs, model, cluster_count ,centroids=None):
    if centroids is None:
       centroids = _generate_kplusplus(contigs,model,cluster_count,DNA.kmer_hash_count)
    cluster_different = True
    while (cluster_different):

        clusters = _expectation(contigs,model,centroids)
#        #Expectations
#        clusters = [set() for _ in xrange(cluster_count)]
#        for contig in contigs:
#            prob = [model.log_probability(contig.signature,centroid) for centroid in centroids]
#            clust_ind = np.argmax(prob)
#            clusters[clust_ind].add(contig)

        new_centroids = _maximization(contigs, model, clusters, centroids.shape)
        #Maximization
#        new_centroids = np.zeros(centroids.shape)
#        for clust_ind ,clust in enumerate(clusters):
#            if not clust:
#                select_as_centroid = np.random.randint(0,len(contigs))
#                new_centroid = model.fit_nonzero_parameters(contigs[select_as_centroid].signature,DNA.kmer_hash_count)
#            else:
#                new_centroid_count = Counter()
#                [new_centroid_count.update(contig.signature) for contig in clust]
#                new_centroid = model.fit_nonzero_parameters(new_centroid_count,DNA.kmer_hash_count)
#            new_centroids[clust_ind,:] = new_centroid
            
        #Keep looping?
        if (new_centroids == centroids).all():
            cluster_different = False
        (new_centroids,centroids) = (centroids,new_centroids)
    return (centroids,clusters)

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
    
def _evaluate_clustering(centroids,clusters):
    cluster_prob = 0
    for i,cluster in enumerate(clusters):
        cluster_prob += sum([model.log_probability(contig.signature,centroids[i]) for contig in cluster])
    return cluster_prob


