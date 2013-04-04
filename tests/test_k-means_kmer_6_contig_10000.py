#!/usr/bin/env python
from nose.tools import  assert_equal
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
from probin.model.composition import multinomial
import numpy as np
import sys

class TestKmeans10000(object):
    FASTA="tests/generated_contigs_10000_test.fna"
    def setUp(self):
        dna.DNA.generate_kmer_hash(6)
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
        contigs_hash = {0:[0,1,2],1:[3,4,5],2:[6,7,8],3:[9,10,11,12,13],4:[14,15,16,17,18],5:[19,20,21,22,23],6:[24,25,26,27,28]}

        for i in xrange(self.cluster_count):
            contig_inds =  contigs_hash[i]
            correct_centroids[i,:] = multinomial.fit_nonzero_parameters([self.contigs[j] for j in contig_inds])


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

    def test_generate_kplusplus_centroids10000(self):
        
        centroids = kmeans._generate_kplusplus(self.contigs,multinomial,self.cluster_count,dna.DNA.kmer_hash_count)
        assert_equal(len(centroids), self.cluster_count)
        assert_equal(len(centroids[0]),dna.DNA.kmer_hash_count )
        assert_equal(np.sum(centroids,axis=1).all(),1)
        
        (clusters,clust_prob,new_centroids) = kmeans.cluster(self.contigs,multinomial,self.cluster_count,max_iter=100,repeat=10)
        
