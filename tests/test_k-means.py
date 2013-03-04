#!/usr/bin/env python
from nose.tools import assert_almost_equal, assert_equal
from probin.dna import DNA
from probin.binning import kmeans
from Bio import SeqIO
import tempfile
import sys
FASTA=""">genome1
GGGGCCCCTTTTTAAAATTATATGCGCGCGCAACACGG
>genome2
ATTATATATGAGAGCGCGCGCGGTGTGTCTCTGCTGC
"""


def test_simple_binning():
    with tempfile.NamedTemporaryFile() as fh:
        fh.write(FASTA)
        fh.seek(0)
        seqs = SeqIO.parse(fh,"fasta")
        seqs = list(seqs)
    contigs = []
    DNA.generate_kmer_hash(4)
    
    print len(DNA.kmer_hash.keys())
    print len(DNA.kmer_hash.values())
    for seq in seqs:
        contigs.append(DNA(seq.id,seq.seq.tostring()))
    assert_equal(True, kmeans.cluster(contigs,2))
    

