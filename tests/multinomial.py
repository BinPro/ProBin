#!/usr/bin/env python
from probin.model.composition import multinomial as ml
import random
import fileinput
import unittest
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from numpy import array

# testing function: log_probability
def test_uniform_one_contig_prob():
    f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
    c = list(SeqIO.parse(f,"fasta"))
    f.close
    signatures = ml.calculate_signatures(4, c)
    k = 4**4
    uniform_prob = [1.0/k]*k
    s=signatures[0]
    log_prob = ml.log_probability(s,uniform_prob)
    assert_almost_equal(log_prob, -1001.20751357)

# testing function: calculate_signatures
def test_signatures_one_contig_basic():
    f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
    c = list(SeqIO.parse(f,"fasta"))
    f.close
    correct_signatures = CORRECT_SIGNATURES_ONE_CONTIG
    calculated_signatures = list(ml.calculate_signatures(4, c)[0])
    assert_equal(calculated_signatures,correct_signatures)

# testing function: fit_parameters
def test_signatures_one_contig_basic():
    f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
    c = list(SeqIO.parse(f,"fasta"))
    f.close
    c_sig = CORRECT_SIGNATURES_ONE_CONTIG
    n = sum(c_sig)
    correct_parameters = [[sig/float(n) for sig in c_sig]]
    calculated_parameters = ml.fit_parameters(4,c)
    assert_equal(calculated_parameters, correct_parameters)

def test_signaturs_large_genome():
    f = fileinput.input("data/8M_genome.fna")
    c= list(SeqIO.parse(f,"fasta"))
    f.close
    calculated_parameters = ml.fit_parameters(4,c)
    assert_equal(len(calculated_parameters[0]), 256)



## "Constants"
CORRECT_SIGNATURES_ONE_CONTIG = [52, 40, 27, 32, 29, 36, 20, 35, 21, 29, 12, 29, 32, 26, 12, 35, 27, 24, 15, 9, 24, 33, 22, 20, 18, 17, 17, 17, 42, 23, 23, 24, 29, 25, 17, 12, 22, 27, 21, 22, 12, 13, 15, 8, 28, 18, 16, 19, 31, 17, 28, 14, 31, 16, 18, 11, 13, 3, 13, 16, 27, 15, 11, 22, 31, 25, 15, 24, 12, 22, 15, 29, 16, 18, 3, 15, 16, 16, 9, 16, 26, 22, 15, 20, 28, 21, 22, 28, 21, 21, 16, 17, 28, 22, 15, 17, 20, 22, 16, 18, 5, 24, 22, 18, 17, 26, 15, 21, 23, 20, 13, 11, 39, 19, 23, 22, 19, 18, 18, 15, 22, 7, 14, 25, 18, 5, 13, 22, 28, 18, 17, 19, 23, 23, 25, 17, 22, 17, 14, 10, 17, 18, 10, 8, 17, 13, 7, 13, 22, 23, 16, 18, 14, 15, 24, 18, 16, 10, 14, 10, 14, 19, 19, 11, 11, 24, 19, 5, 20, 7, 9, 13, 21, 14, 11, 5, 39, 18, 28, 13, 13, 20, 20, 14, 19, 12, 10, 7, 23, 10, 10, 8, 40, 37, 32, 30, 10, 18, 9, 31, 24, 28, 19, 27, 25, 16, 14, 16, 25, 19, 15, 15, 9, 22, 15, 16, 23, 16, 22, 15, 17, 15, 16, 7, 20, 22, 11, 12, 12, 4, 9, 5, 14, 13, 10, 9, 26, 15, 8, 16, 30, 14, 19, 22, 11, 8, 20, 15, 11, 8, 9, 17, 17, 24, 11, 15]
