#!/usr/bin/env python
import os
import sys
import fileinput
from nose.tools import assert_equal, assert_true
import numpy as np
from Bio import SeqIO

from probin import dna
from probin.binning import kmeans
from probin.binning import statistics
from probin.model.composition import multinomial

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))

class TestBinningStats(object):
    FASTA="tests/generated_contigs_10000_test.fna"
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
        self.cluster_count = 3
        fh = fileinput.input(os.path.join(data_path,"generated_contigs_10000_test.fna"))
        seqs = SeqIO.parse(fh,"fasta")
        seqs = list(seqs)
        self.contigs = []

        for seq in seqs:
            self.contigs.append(dna.DNA(seq.description,seq.seq.tostring(),calc_sign=True))
        (self.clusters,self.clust_prob,self.centroids) = kmeans.cluster(self.contigs, multinomial.log_probability, multinomial.fit_nonzero_parameters, self.cluster_count)

    def tearDown(self):
        reload(dna)
        self.contigs = []
            
    def test_recall(self):
        assert_equal(True,True)
