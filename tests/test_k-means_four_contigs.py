#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal
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

class TestKmeans(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
        self.cluster_count = 2
        fh = fileinput.input(os.path.join(data_path,"generated_contigs_10000_test.fna"))
        seqs = SeqIO.parse(fh,"fasta")
        seqs = list(seqs)
        self.contigs = []
        self.rs = np.random.RandomState(seed=1)
        
        for seq in seqs:
            if  seq.id.startswith("Ehrlichia_canis_Jake_uid58071_1001") or \
                seq.id.startswith("Ehrlichia_canis_Jake_uid58071_1002") or \
                seq.id.startswith("Ehrlichia_ruminantium_Welgevonden_uid58243_1033") or \
                seq.id.startswith("Ehrlichia_ruminantium_Welgevonden_uid58243_1034"):

                self.contigs.append(dna.DNA(seq.id,seq.seq.tostring()))
            
        for contig in self.contigs:
            contig.calculate_signature()
            

    def tearDown(self):
        reload(dna)
        self.contigs = []

    def test_contigs_created(self):
        c_ids = set(["Ehrlichia_canis_Jake_uid58071_1001", "Ehrlichia_canis_Jake_uid58071_1002","Ehrlichia_ruminantium_Welgevonden_uid58243_1033", "Ehrlichia_ruminantium_Welgevonden_uid58243_1034"])
        for contig in self.contigs:
            assert_true(contig.id in c_ids)
            
    def test_generate_centroids(self):
        centroids = kmeans._generate_centroids(self.cluster_count,5,self.rs)
        assert_equal(len(centroids), self.cluster_count)
        assert_equal(len(centroids[0]),5 )
        assert_equal(np.sum(centroids,axis=1).all(),1)

    def test_generate_kplusplus_centroids(self):
        centroids = kmeans._generate_kplusplus(self.contigs,multinomial.log_probability, multinomial.fit_nonzero_parameters,self.cluster_count,dna.DNA.kmer_hash_count,self.rs)
        assert_equal(len(centroids), self.cluster_count)
        assert_equal(len(centroids[0]),dna.DNA.kmer_hash_count )
        assert_equal(np.sum(centroids,axis=1).all(),1)
        
        correct_centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        correct_centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0], self.contigs[1]])
        correct_centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2], self.contigs[3]])
        correct_clusters = kmeans._expectation(self.contigs,multinomial.log_probability,correct_centroids)
        correct_clust_prob = kmeans._evaluate_clustering(multinomial.log_probability, correct_clusters,correct_centroids)
        
        (clusters, clust_prob,new_centroids) = kmeans.cluster(self.contigs,multinomial.log_probability, multinomial.fit_nonzero_parameters,self.cluster_count,centroids)
        print clust_prob
        assert_almost_equal(0, min(np.abs(clust_prob - np.array([-1659.9510320847476, -1652.322663414292, -1658.28785337, -1665.52431153]))))
        
    def test_cluster_perfect_center(self):
        centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0],self.contigs[1]])
        centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2],self.contigs[3]])
        correct_clusters = kmeans._expectation(self.contigs,multinomial.log_probability,centroids)        

        (clusters, clust_prob,new_centroids) = kmeans.cluster(self.contigs,multinomial.log_probability, multinomial.fit_nonzero_parameters,self.cluster_count,centroids)
        
        assert_equal(kmeans._evaluate_clustering(multinomial.log_probability, correct_clusters, centroids),clust_prob)

    def test_cluster_semi_center(self):
        centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0]])
        centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2]])
        (clusters,clust_prob,new_centroids) = kmeans.cluster(self.contigs,multinomial.log_probability, multinomial.fit_nonzero_parameters,2,centroids)

        correct_centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        correct_centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0], self.contigs[1]])
        correct_centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2], self.contigs[3]])
        correct_clusters = kmeans._expectation(self.contigs,multinomial.log_probability,correct_centroids)        
        
        assert_equal(kmeans._evaluate_clustering(multinomial.log_probability, correct_clusters, correct_centroids),clust_prob)
        
        
