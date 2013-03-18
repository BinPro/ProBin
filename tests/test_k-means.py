#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal, assert_true
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
import tempfile
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
        centroids = kmeans._generate_centroids(2,5)
        assert_equal(len(centroids), 2)
        assert_equal(len(centroids[0]),5 )
        assert_equal(np.sum(centroids,axis=1).all(),1)
        
    def test_cluster_perfect_center(self):
        centroids = np.zeros((2,dna.DNA.kmer_hash_count))
        centroids[0,:] = multinomial.fit_nonzero_parameters(self.contigs[0].signature + self.contigs[1].signature,dna.DNA.kmer_hash_count)
        centroids[1,:] = multinomial.fit_nonzero_parameters(self.contigs[2].signature + self.contigs[3].signature,dna.DNA.kmer_hash_count)
        new_centroids = kmeans.cluster(self.contigs,2,multinomial,centroids)
        assert_equal((centroids==new_centroids).all(),True)

    def test_cluster_semi_center(self):
        centroids = np.zeros((2,dna.DNA.kmer_hash_count))
        centroids[0,:] = multinomial.fit_nonzero_parameters(self.contigs[0].signature,dna.DNA.kmer_hash_count)
        centroids[1,:] = multinomial.fit_nonzero_parameters(self.contigs[2].signature,dna.DNA.kmer_hash_count)
        new_centroids = kmeans.cluster(self.contigs,2,multinomial,centroids)

        correct_centroids = np.zeros((2,dna.DNA.kmer_hash_count))
        correct_centroids[0,:] = multinomial.fit_nonzero_parameters(self.contigs[0].signature + self.contigs[1].signature,dna.DNA.kmer_hash_count)
        correct_centroids[1,:] = multinomial.fit_nonzero_parameters(self.contigs[2].signature + self.contigs[3].signature,dna.DNA.kmer_hash_count)        

        assert_equal((correct_centroids==new_centroids).all(),True)
        
        
        