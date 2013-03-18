#!/usr/bin/env python
from probin.model.composition import multinomial as ml
from probin import dna
import fileinput
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from collections import Counter
import numpy as np

class TestMultinomial(object):
    ## "Constants"
    CORRECT_SIGNATURES_ONE_CONTIG = Counter({0: 73, 1: 73, 27: 64, 73: 62, 75: 60, 15: 59, 21: 59, 58: 59, 7: 58, 12: 56, 74: 56, 77: 56, 102: 56, 3: 55, 13: 53, 55: 53, 45: 52, 63: 52, 4: 51, 31: 51, 2: 49, 11: 49, 16: 49, 20: 49, 48: 49, 28: 48, 41: 48, 78: 48, 83: 47, 122: 47, 126: 47, 49: 46, 60: 46, 76: 46, 91: 46, 110: 46, 111: 46, 9: 45, 29: 45, 46: 44, 87: 44, 97: 44, 30: 43, 121: 43, 93: 42, 99: 42, 131: 41, 36: 40, 71: 40, 8: 39, 39: 39, 89: 39, 128: 39, 6: 38, 22: 38, 34: 38, 23: 37, 35: 37, 72: 37, 94: 37, 101: 37, 5: 36, 92: 36, 26: 35, 64: 35, 79: 35, 84: 35, 85: 35, 116: 35, 18: 34, 24: 34, 82: 34, 44: 33, 81: 33, 100: 33, 32: 32, 86: 32, 90: 32, 25: 31, 43: 31, 80: 31, 112: 31, 47: 29, 57: 29, 104: 29, 38: 28, 50: 28, 54: 28, 66: 28, 69: 27, 98: 27, 119: 27, 123: 27, 132: 27, 70: 26, 113: 26, 115: 26, 62: 25, 96: 25, 17: 24, 56: 24, 59: 24, 105: 24, 117: 24, 88: 23, 129: 23, 19: 22, 33: 22, 95: 22, 108: 22, 53: 21, 67: 21, 68: 21, 114: 21, 118: 21, 133: 21, 10: 20, 51: 20, 125: 19, 40: 18, 42: 18, 37: 17, 103: 17, 14: 16, 106: 15, 130: 15, 65: 14, 107: 13, 109: 13, 61: 12, 127: 9, 135: 9, 134: 8, 124: 7, 120: 5, 52: 3})
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    def tearDown(self):
        reload(dna)
    # testing function: log_probability
    def test_uniform_one_contig_prob(self):
        f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
        c = list(SeqIO.parse(f,"fasta"))
        f.close()
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()
        s = dna_c.signature
        k = 4**4
        uniform_prob = {}
        for i,cnt in s.items():
            uniform_prob[i] = 1./k
        log_prob = ml.log_probability(s,uniform_prob)
        print log_prob
        assert_almost_equal(log_prob, -3791.05738056)
    
    # testing function: calculate_signatures
    def test_signatures_one_contig_basic(self):
        f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
        c = list(SeqIO.parse(f,"fasta"))
        f.close()
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()
        calculated_signature = dna_c.signature
        correct_signature = self.CORRECT_SIGNATURES_ONE_CONTIG
        assert_equal(calculated_signature,correct_signature)
    
    # testing function: fit_parameters
    def test_parameters_one_contig_basic(self):
        f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
        c = list(SeqIO.parse(f,"fasta"))
        f.close()
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()
        c_sig = self.CORRECT_SIGNATURES_ONE_CONTIG
        n = sum(c_sig.values())
        correct_parameters = {}
        for i,v in c_sig.items():
            correct_parameters[i] = v/float(n)
        calculated_parameters = ml.fit_parameters(dna_c.signature)
        assert_equal(calculated_parameters, correct_parameters)
    
    def test_signaturs_large_genome(self):
        f = fileinput.input("data/8M_genome.fna")
        c= list(SeqIO.parse(f,"fasta"))
        f.close
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()
        calculated_parameters = ml.fit_parameters(dna_c.signature)
        assert_equal(len(calculated_parameters), 136)


    def test_fit_nonzero_parameters(self):
        c = Counter([1,2,2,3,3,3])
        distribution = ml.fit_nonzero_parameters(c,6)
        pseudo = np.ones(6)
        true_dist = (pseudo + np.array([0,1,2,3,0,0]))/float(12)
        assert_equal((true_dist==distribution).all(), True)

    def test_log_probability(self):
        c = Counter({0:10,1: 4,2:1,3:8,4:20})
        log_p_val = ml.log_probability(c,[0.2,0.1,0.05,0.2,0.45])
        real_p = np.exp(log_p_val)
        assert_almost_equal(real_p, 0.001074701)

        c = Counter({0:1, 1:0})
        log_p_val = ml.log_probability(c,[0.5,0.5])
        real_p = np.exp(log_p_val)
        assert_almost_equal(real_p, 0.5)

        c = Counter({0:0, 1:1})
        log_p_val = ml.log_probability(c,[0.5,0.5])
        real_p = np.exp(log_p_val)
        assert_almost_equal(real_p, 0.5)
