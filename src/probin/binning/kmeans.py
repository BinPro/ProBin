#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from probin.dna import DNA
import numpy as np
from numpy.random import randint
from collections import Counter

def cluster(contigs, cluster_count,model):
    centroids = _generate_centroids(contigs,cluster_count,model)
    
    clusters = [set() for i in xrange(cluster_count)]
    
    keep_clustering = True
    counter = 0
    while keep_clustering:
        keep_clustering=False
        clusters_new = _expectation(contigs,centroids,model)
        centroids = _maximization(clusters,model)
        number_same = 0
        for cluster in clusters_new:
            if cluster not in clusters:
                keep_clustering = True
            else:
                number_same += 1
            if cluster == set():
                print "empty"
        (clusters, clusters_new) = (clusters_new,clusters)
        counter += 1
        print number_same
        print counter
        print "-"*70
    return clusters
    
def _generate_centroids(contigs,cluster_count ,model):
    centroids = np.zeros((cluster_count,DNA.kmer_hash_count))
    centroid_contigs = set()
    first_centroid_contig = randint(low=0,high=cluster_count)
    centroid_contigs.add(first_centroid_contig)
    print first_centroid_contig
    centroids[0,:] = model.fit_nonzero_parameters(contigs[first_centroid_contig].signature,DNA.kmer_hash_count)
    for ind_next_centroid in xrange(1,cluster_count):
        dist = np.zeros(len(contigs))
        for ind_contig,contig in enumerate(contigs):
            for ind_done_centroid in xrange(ind_next_centroid):
                dist[ind_contig] = min(dist[ind_contig],model.log_probability(contig.signature,centroids[ind_done_centroid]))
        dist = np.power(dist,2)
        dist /= np.sum(dist)
        while True:
            next_centroid_contig = np.argmax(np.random.random(dist.shape) * dist)
            if next_centroid_contig not in centroid_contigs:
                centroid_contigs.add(next_centroid_contig)
                break
            print "had to find another"
        print next_centroid_contig
        centroids[ind_next_centroid,:] = model.fit_nonzero_parameters(contigs[next_centroid_contig].signature,DNA.kmer_hash_count)
    print "done centroids"
    return centroids
    
def _expectation(contigs,centroids,model):
    clusters_new = [set() for i in xrange(len(centroids))]
    for i,contig in enumerate(contigs):
        max_probability = None
        max_clusters = []        
        for centroid_index,centroid in enumerate(centroids):
            centroid_probability = model.log_probability(contig.signature,centroid)

            if centroid_probability >= max_probability:
                if centroid_probability > max_probability:
                    max_clusters = []
                max_probability = centroid_probability
                max_clusters.append(centroid_index)
        clusters_new[max_clusters[0]].add(contig)
    return clusters_new
    
def _maximization(clusters,model):
    centroids_new = np.zeros((len(clusters),DNA.kmer_hash_count))
    for i,cluster in enumerate(clusters):
        if not len(cluster):
            c = np.random.rand(DNA.kmer_hash_count)
            c /= np.sum(c)
            centroids_new[i,:] = c
        else:        
            clust_sig = Counter()
            [clust_sig.update(contig.signature) for contig in cluster]
            centroids_new[i,:] = model.fit_nonzero_parameters(clust_sig,DNA.kmer_hash_count)
    return centroids_new

if __name__=="__main__":
    from Bio import SeqIO
    from probin.model.composition import multinomial
    FASTA="../../../tests/generated_contigs_test.fna"
    DNA.generate_kmer_hash(4)
    cluster_count = 2
    with open(FASTA,"r") as fh:
        seqs = SeqIO.parse(fh,"fasta")
        seqs = list(seqs)
    contigs = []

    for seq in seqs:
        contigs.append(DNA(seq.id,seq.seq.tostring()))
    for contig in contigs:
        contig.calculate_signature()
    clusters = cluster(contigs,7,multinomial)

