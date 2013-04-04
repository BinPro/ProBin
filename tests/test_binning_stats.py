#!/usr/bin/env python
import os
import sys
from nose.tools import assert_equal, assert_true
import numpy as np
from Bio import SeqIO

from probin import dna
from probin.binning import kmeans
from probin.binning import statistics
from probin.model.composition import multinomial


class TestBinningStats(object):
    FASTA="tests/generated_contigs_10000_test.fna"
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
        self.cluster_count = 3
        with open(self.FASTA,"r") as fh:
            seqs = SeqIO.parse(fh,"fasta")
            seqs = list(seqs)
        self.contigs = []

        for seq in seqs:
            self.contigs.append(dna.DNA(seq.description,seq.seq.tostring(),calc_sign=True))
        for contig in self.contigs:
            (contig.phylo_family, contig.phylo_genus, contig.phylo_species) = contig.id.split(" ",1)[1].split("|")
        (self.clusters,self.clust_prob,self.centroids) = kmeans.cluster(self.contigs, multinomial, self.cluster_count)

    def tearDown(self):
        reload(dna)
        self.contigs = []
            
    def test_recall(self):
        recall = statistics.recall(self.clusters,self.contigs)
        print>>sys.stderr, recall