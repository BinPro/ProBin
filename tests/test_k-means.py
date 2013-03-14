#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
import tempfile
from probin.model.composition import multinomial
import numpy as np
import sys

class TestKmeans(object):
    FASTA=""">genome1
GGGGCCCCTTTTTAAAATTATATGCGCGCGACAACACTG
>genome2
ATTATATATGAGAGATACGCGCGCTGTGTCTCTGCTGC
>genome3
GGGGCCCCTTTTTAAAATTATATGCGCGCGCAACACGG
>genome4
ATTATATATGAGAGCGCGCGCGGTGTGTCTCTGCTGC
"""
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
        self.cluster_count = 2
        with tempfile.NamedTemporaryFile() as fh:
            fh.write(self.FASTA)
            fh.seek(0)
            seqs = SeqIO.parse(fh,"fasta")
            seqs = list(seqs)
        self.contigs = []

        for seq in seqs:
            self.contigs.append(dna.DNA(seq.id,seq.seq.tostring()))
        for contig in self.contigs:
            contig.calculate_signature()
            

    def tearDown(self):
        reload(dna)
        self.contigs = []
        
    def test_simple_binning(self):
        clusters = kmeans.cluster(self.contigs,self.cluster_count,multinomial)
        s = set()
        [set.add(x) for x in clusters]
        assert_equal(len(self.contigs),len(s)) 
        
    def test_expectation(self):
        centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        centroids[0,:] = multinomial.fit_nonzero_parameters(self.contigs[0].signature,dna.DNA.kmer_hash_count)
        centroids[1,:] = multinomial.fit_nonzero_parameters(self.contigs[2].signature,dna.DNA.kmer_hash_count)
        kmeans._expectation(self.contigs,centroids,multinomial)
        
    
    
    def test_maximization(self):
        pass
    