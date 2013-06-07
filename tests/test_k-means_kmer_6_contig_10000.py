#!/usr/bin/env python
from nose.tools import  assert_equal
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
from probin.model.composition import multinomial
import numpy as np
import sys
import os
import fileinput

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))

class TestKmeans10000(object):
    def setUp(self):
        reload(dna)
        reload(kmeans)
        dna.DNA.generate_kmer_hash(6)
        print >> sys.stderr, dna.DNA.kmer_hash_count
        self.cluster_count = 7
        fh = fileinput.input(os.path.join(data_path,"generated_contigs_10000_test.fna"))
        seqs = SeqIO.parse(fh,"fasta")
        seqs = list(seqs)
        self.contigs = []

        for seq in seqs:
            self.contigs.append(dna.DNA(seq.id,seq.seq.tostring()))
            
        for contig in self.contigs:
            contig.calculate_signature()

        correct_centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        contigs_hash = {0:[0,1,2],1:[3,4,5],2:[6,7,8],3:[9,10,11,12,13],4:[14,15,16,17,18],5:[19,20,21,22,23],6:[24,25,26,27,28]}

        for i in xrange(self.cluster_count):
            contig_inds =  contigs_hash[i]
            correct_centroids[i,:] = multinomial.fit_nonzero_parameters([self.contigs[j] for j in contig_inds])
        
        self.correct_centroids = correct_centroids
        self.correct_clusters = kmeans._expectation(self.contigs,multinomial.log_probabilities,self.correct_centroids)        
        self.rs = np.random.RandomState(seed=1)
	self.params = {"contigs":self.contigs}
	self.max_iter=100
	self.run = 3
	self.epsilon = 0.001
	self.verbose = False
    def tearDown(self):
        reload(dna)
        reload(kmeans)
        self.contigs = []

    def test_generate_kplusplus_centroids10000(self):
        
        centroids = kmeans._generate_kplusplus(self.contigs,multinomial.log_probabilities,multinomial.fit_nonzero_parameters,self.cluster_count,dna.DNA.kmer_hash_count,self.rs)
        print>>sys.stderr, centroids.shape
        print>>sys.stderr, len(self.contigs)
        print>>sys.stderr, self.cluster_count
	self.params["centroids"] = centroids
        (clusters,clust_prob,new_centroids) = kmeans._clustering(self.cluster_count,self.max_iter, self.run, self.epsilon, self.verbose, multinomial.log_probabilities,multinomial.fit_nonzero_parameters,**self.params)
        assert_equal(len(centroids), self.cluster_count)
        assert_equal(len(centroids[0]),dna.DNA.kmer_hash_count )
        assert_equal(np.sum(centroids,axis=1).all(),1)
        
        
