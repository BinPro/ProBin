#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal
from probin import dna
from probin.binning import kmeans
from Bio import SeqIO
import tempfile
import sys


class KmeansTest(object):
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
        
        print len(DNA.kmer_hash.keys())
        print max(DNA.kmer_hash.values())
        for seq in seqs:
            contigs.append(DNA(seq.id,seq.seq.tostring()))
        assert_equal(True, kmeans.cluster(contigs,2))
    

