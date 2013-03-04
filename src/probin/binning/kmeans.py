#!/usr/bin/env python
"""Cluster DNA based on k-means algorithm with fixed number of clusters"""
from sklearn.cluster import KMeans
import sys
def cluster(contigs, cluster_count):
    estimator = KMeans(cluster_count)
    print>>sys.stderr, contigs[0].kmer_hash_count
    estimator.fit([contig.signature for contig in contigs])
    


