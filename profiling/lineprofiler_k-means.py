#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:11:45 2013

@author: binni
"""


from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
from probin.model.composition import multinomial
import sys

FASTA="../tests/generated_contigs_10000_test.fna"
def setUp():
    dna.DNA.generate_kmer_hash(6)
    print>>sys.stderr,dna.DNA.kmer_hash_count
    cluster_count = 7
    with open(FASTA,"r") as fh:
        seqs = SeqIO.parse(fh,"fasta")
        seqs = list(seqs)
    contigs = []

    for seq in seqs:
        contigs.append(dna.DNA(seq.id,seq.seq.tostring()))
        
    for contig in contigs:
        contig.calculate_signature()

    return (contigs,cluster_count)


def test_generate_kplusplus_centroids():
    print>>sys.stderr,"1, kplusplus_cenroids"
    (contigs,cluster_count) = setUp()
    print>>sys.stderr,"2,kplusplus_cenroids" 
    
    (clust_prob,new_centroids,clusters) = kmeans.cluster(contigs,multinomial,cluster_count,max_iter=3,repeat=2)
    
    print>>sys.stderr,"3,kplusplus_cenroids"        
    for c in clusters:
        print>>sys.stderr,c
    print "-"*70

if __name__=="__main__":
    test_generate_kplusplus_centroids()