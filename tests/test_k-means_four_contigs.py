#!/usr/bin/env python
from nose.tools import assert_equal, assert_true
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
from probin.model.composition import multinomial
import numpy as np
import sys

class TestKmeans(object):
    FASTA="tests/generated_contigs_test.fna"
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
        self.cluster_count = 2
        with open(self.FASTA,"r") as fh:
            seqs = SeqIO.parse(fh,"fasta")
            seqs = list(seqs)
        self.contigs = []

        for seq in seqs:
            if  seq.id.startswith("Ehrlichia_canis_Jake_uid58071_1000") or \
                seq.id.startswith("Ehrlichia_canis_Jake_uid58071_1001") or \
                seq.id.startswith("Ehrlichia_ruminantium_Welgevonden_uid58243_1003") or \
                seq.id.startswith("Ehrlichia_ruminantium_Welgevonden_uid58243_1004"):

                self.contigs.append(dna.DNA(seq.id,seq.seq.tostring()))
            
        for contig in self.contigs:
            contig.calculate_signature()
            

    def tearDown(self):
        reload(dna)
        self.contigs = []

    def test_contigs_created(self):
        c_ids = set(["Ehrlichia_canis_Jake_uid58071_1000", "Ehrlichia_canis_Jake_uid58071_1001","Ehrlichia_ruminantium_Welgevonden_uid58243_1003", "Ehrlichia_ruminantium_Welgevonden_uid58243_1004"])
        for contig in self.contigs:
            assert_true(contig.id in c_ids)
            
    def test_generate_centroids(self):
        centroids = kmeans._generate_centroids(self.cluster_count,5)
        assert_equal(len(centroids), self.cluster_count)
        assert_equal(len(centroids[0]),5 )
        assert_equal(np.sum(centroids,axis=1).all(),1)

    def test_generate_kplusplus_centroids(self):
        centroids = kmeans._generate_kplusplus(self.contigs,multinomial,2,dna.DNA.kmer_hash_count)
        assert_equal(len(centroids), self.cluster_count)
        assert_equal(len(centroids[0]),dna.DNA.kmer_hash_count )
        assert_equal(np.sum(centroids,axis=1).all(),1)
        
        correct_centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        correct_centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0], self.contigs[1]])
        correct_centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2], self.contigs[3]])
        correct_clusters = kmeans._expectation(self.contigs,multinomial,correct_centroids)        
        
        (clusters, clust_prob,new_centroids) = kmeans.cluster(self.contigs,multinomial,2,centroids)
        
        assert_equal(kmeans._evaluate_clustering(multinomial, correct_clusters, correct_centroids),clust_prob)

                        
        
    def test_cluster_perfect_center(self):
        centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0],self.contigs[1]])
        centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2],self.contigs[3]])
        correct_clusters = kmeans._expectation(self.contigs,multinomial,centroids)        

        (clusters, clust_prob,new_centroids) = kmeans.cluster(self.contigs,multinomial,2,centroids)
        
        assert_equal(kmeans._evaluate_clustering(multinomial, correct_clusters, centroids),clust_prob)

    def test_cluster_semi_center(self):
        centroids = np.zeros((2,dna.DNA.kmer_hash_count))
        centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0]])
        centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2]])
        (clusters,clust_prob,new_centroids) = kmeans.cluster(self.contigs,multinomial,2,centroids)

        correct_centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        correct_centroids[0,:] = multinomial.fit_nonzero_parameters([self.contigs[0], self.contigs[1]])
        correct_centroids[1,:] = multinomial.fit_nonzero_parameters([self.contigs[2], self.contigs[3]])
        correct_clusters = kmeans._expectation(self.contigs,multinomial,correct_centroids)        
        
        assert_equal(kmeans._evaluate_clustering(multinomial, correct_clusters, correct_centroids),clust_prob)
        
        
