#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
import tempfile
import sys


class TestKmeans(object):
    FASTA=""">genome1
GGGGCCCCTTTTTAAAATTATATGCGCGCGCAACACGG
>genome2
ATTATATATGAGAGCGCGCGCGGTGTGTCTCTGCTGC
"""
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    
    def tearDown(self):
        reload(dna)

    def test_simple_binning(self):
        with tempfile.NamedTemporaryFile() as fh:
            fh.write(self.FASTA)
            fh.seek(0)
            seqs = SeqIO.parse(fh,"fasta")
            seqs = list(seqs)
        contigs = []
        
        print len(dna.DNA.kmer_hash.keys())
        print max(dna.DNA.kmer_hash.values())
        for seq in seqs:
            contigs.append(dna.DNA(seq.id,seq.seq.tostring()))
        assert_equal(True,True)
    

