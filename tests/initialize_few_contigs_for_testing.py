# -*- coding: utf-8 -*-
"""
Created on Thu May  2 17:05:21 2013

@author: binni

Initialize few contigs for testing
"""
import sys
import os

from probin.dna import DNA

from probin.model.composition import multinomial
from probin.binning import em,kmeans
from probin.binning import statistics

DNA.generate_kmer_hash(4)
cluster_count=3

contigs = [DNA(id=1, seq="AGAGAGAGATATATAT",calc_sign=True),
           DNA(id=2, seq="AGAGAGACATATATAT",calc_sign=True),
           DNA(id=3, seq="ACACACACTGTGTGTG",calc_sign=True),
           DNA(id=4, seq="ACACACATTGTGTGTG",calc_sign=True),
           DNA(id=5, seq="GACGACGACTACTACT",calc_sign=True),
           DNA(id=6, seq="GACGACGAGTACTACT",calc_sign=True)]
phylos = ["A1|B1|C1","A1|B1|C1","A2|B2|C2","A2|B2|C2","A3|B3|C3","A3|B3|C3"]
for (c,p) in zip(contigs,phylos):
    c.phylo = p

print kmeans.cluster(contigs,multinomial,cluster_count,centroids=None,max_iter=100, repeat=10)
print em.cluster(contigs,multinomial,cluster_count,centroids=None,max_iter=100, repeat=10,epsilon=0.01)