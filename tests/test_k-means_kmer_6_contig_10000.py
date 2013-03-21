#!/usr/bin/env python
from nose.tools import  assert_equal
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
from probin.model.composition import multinomial
import numpy as np
import sys

class TestKmeans(object):
    FASTA="tests/generated_contigs_10000_test.fna"
    def setUp(self):
        dna.DNA.generate_kmer_hash(6)
        print>>sys.stderr,dna.DNA.kmer_hash_count
        self.cluster_count = 7
        with open(self.FASTA,"r") as fh:
            seqs = SeqIO.parse(fh,"fasta")
            seqs = list(seqs)
        self.contigs = []

        for seq in seqs:
            self.contigs.append(dna.DNA(seq.id,seq.seq.tostring()))
            
        for contig in self.contigs:
            contig.calculate_signature()

        correct_centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
        correct_centroids[0,:] = multinomial.fit_nonzero_parameters(self.contigs[0].signature + self.contigs[1].signature + self.contigs[2].signature,dna.DNA.kmer_hash_count)
        correct_centroids[1,:] = multinomial.fit_nonzero_parameters(self.contigs[3].signature + self.contigs[4].signature + self.contigs[5].signature,dna.DNA.kmer_hash_count)
        correct_centroids[2,:] = multinomial.fit_nonzero_parameters(self.contigs[6].signature + self.contigs[7].signature + self.contigs[8].signature,dna.DNA.kmer_hash_count)
        correct_centroids[3,:] = multinomial.fit_nonzero_parameters(self.contigs[9].signature + self.contigs[10].signature + self.contigs[11].signature + self.contigs[12].signature + self.contigs[13].signature,dna.DNA.kmer_hash_count)
        correct_centroids[4,:] = multinomial.fit_nonzero_parameters(self.contigs[14].signature + self.contigs[15].signature + self.contigs[16].signature + self.contigs[17].signature + self.contigs[18].signature,dna.DNA.kmer_hash_count)
        correct_centroids[5,:] = multinomial.fit_nonzero_parameters(self.contigs[19].signature + self.contigs[20].signature + self.contigs[21].signature + self.contigs[22].signature + self.contigs[23].signature,dna.DNA.kmer_hash_count)
        correct_centroids[6,:] = multinomial.fit_nonzero_parameters(self.contigs[24].signature + self.contigs[25].signature + self.contigs[26].signature + self.contigs[27].signature + self.contigs[28].signature,dna.DNA.kmer_hash_count)
#        correct_centroids[0,:] = multinomial.fit_nonzero_parameters(\
#                                     self.contigs[0].signature + self.contigs[1].signature + self.contigs[2].signature +\
#                                     self.contigs[3].signature + self.contigs[4].signature + self.contigs[5].signature +\
#                                     self.contigs[6].signature + self.contigs[7].signature + self.contigs[8].signature \
#                                     ,dna.DNA.kmer_hash_count)
#        correct_centroids[1,:] = multinomial.fit_nonzero_parameters(\
#                                     self.contigs[9].signature + self.contigs[10].signature + self.contigs[11].signature +\
#                                     self.contigs[12].signature + self.contigs[13].signature + self.contigs[14].signature +\
#                                     self.contigs[15].signature + self.contigs[16].signature + self.contigs[17].signature +\
#                                     self.contigs[18].signature ,dna.DNA.kmer_hash_count)
#        correct_centroids[2,:] = multinomial.fit_nonzero_parameters(\
#                                     self.contigs[19].signature + self.contigs[20].signature + self.contigs[21].signature +\
#                                     self.contigs[22].signature + self.contigs[23].signature + self.contigs[24].signature +\
#                                     self.contigs[25].signature + self.contigs[26].signature + self.contigs[27].signature +\
#                                     self.contigs[28].signature,dna.DNA.kmer_hash_count)

        self.correct_centroids = correct_centroids
        self.correct_clusters = kmeans._expectation(self.contigs,multinomial,self.correct_centroids)        
            

    def tearDown(self):
        reload(dna)
        self.contigs = []

    def test_generate_kplusplus_centroids(self):
        print>>sys.stderr,"1, kplusplus_cenroids"        
        centroids = kmeans._generate_kplusplus(self.contigs,multinomial,self.cluster_count,dna.DNA.kmer_hash_count)
        print>>sys.stderr,"2,kplusplus_cenroids" 
        assert_equal(len(centroids), self.cluster_count)
        assert_equal(len(centroids[0]),dna.DNA.kmer_hash_count )
        assert_equal(np.sum(centroids,axis=1).all(),1)
        
        (clust_prob,new_centroids,clusters) = kmeans.cluster(self.contigs,multinomial,self.cluster_count,max_iter=5,repeat=2)
        
        print>>sys.stderr,"3,kplusplus_cenroids"        
        for c in clusters:
            print>>sys.stderr,c
        print "-"*70
        assert_equal(kmeans._evaluate_clustering(self.correct_centroids,self.correct_clusters,multinomial),clust_prob)

                        
        
#    def test_cluster_perfect_center(self):
#        print>>sys.stderr,"1,perfect_center"
#
#        correct_clusters = kmeans._expectation(self.contigs,multinomial,self.correct_centroids)
#        print>>sys.stderr,"2,perfect_center"
#        (clust_prob,new_centroids,clusters) = kmeans.cluster(self.contigs,multinomial,self.cluster_count,self.correct_centroids,max_iter=5,repeat=2)
#        print>>sys.stderr,"3,perfect_center"
#        for c in clusters:
#            print>>sys.stderr,c
#        print "-"*70
#        assert_equal(kmeans._evaluate_clustering(self.correct_centroids,correct_clusters,multinomial),clust_prob)
#
#    def test_cluster_semi_center(self):
#        print>>sys.stderr,"1,semi_cenroids"        
#
#        correct_clusters = kmeans._expectation(self.contigs,multinomial,self.correct_centroids)
#    
#        semi_correct_centroids = np.zeros((self.cluster_count,dna.DNA.kmer_hash_count))
#        
##        semi_correct_centroids[0,:] = multinomial.fit_nonzero_parameters(self.contigs[0].signature, dna.DNA.kmer_hash_count)
##        semi_correct_centroids[1,:] = multinomial.fit_nonzero_parameters(self.contigs[9].signature, dna.DNA.kmer_hash_count)
##        semi_correct_centroids[2,:] = multinomial.fit_nonzero_parameters(self.contigs[19].signature, dna.DNA.kmer_hash_count)
#
#        semi_correct_centroids[0,:] = multinomial.fit_nonzero_parameters(self.contigs[0].signature, dna.DNA.kmer_hash_count)
#        semi_correct_centroids[1,:] = multinomial.fit_nonzero_parameters(self.contigs[3].signature,dna.DNA.kmer_hash_count)
#        semi_correct_centroids[2,:] = multinomial.fit_nonzero_parameters(self.contigs[6].signature,dna.DNA.kmer_hash_count)
#        semi_correct_centroids[3,:] = multinomial.fit_nonzero_parameters(self.contigs[9].signature,dna.DNA.kmer_hash_count)
#        semi_correct_centroids[4,:] = multinomial.fit_nonzero_parameters(self.contigs[14].signature,dna.DNA.kmer_hash_count)
#        semi_correct_centroids[5,:] = multinomial.fit_nonzero_parameters(self.contigs[19].signature,dna.DNA.kmer_hash_count)
#        semi_correct_centroids[6,:] = multinomial.fit_nonzero_parameters(self.contigs[24].signature,dna.DNA.kmer_hash_count)
#
#        (clust_prob,new_centroids,clusters) = kmeans.cluster(self.contigs,multinomial,self.cluster_count,semi_correct_centroids,max_iter=5,repeat=2)
#        print>>sys.stderr,"3,semi_cenroids"        
#        for c in clusters:
#            print>>sys.stderr,c
#        print "-"*70
#        assert_equal(kmeans._evaluate_clustering(self.correct_centroids,correct_clusters,multinomial),clust_prob)
#        
#        