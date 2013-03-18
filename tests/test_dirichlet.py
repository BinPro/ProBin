#!/usr/bin/env python
from probin.model.composition import dirichlet as model
from probin import dna
import fileinput
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from collections import Counter
import numpy as np

class TestDirichlet(object):
    ## "Constants"
    CORRECT_SIGNATURES_ONE_CONTIG = Counter({0: 73, 1: 73, 27: 64, 73: 62, 75: 60, 15: 59, 21: 59, 58: 59, 7: 58, 12: 56, 74: 56, 77: 56, 102: 56, 3: 55, 13: 53, 55: 53, 45: 52, 63: 52, 4: 51, 31: 51, 2: 49, 11: 49, 16: 49, 20: 49, 48: 49, 28: 48, 41: 48, 78: 48, 83: 47, 122: 47, 126: 47, 49: 46, 60: 46, 76: 46, 91: 46, 110: 46, 111: 46, 9: 45, 29: 45, 46: 44, 87: 44, 97: 44, 30: 43, 121: 43, 93: 42, 99: 42, 131: 41, 36: 40, 71: 40, 8: 39, 39: 39, 89: 39, 128: 39, 6: 38, 22: 38, 34: 38, 23: 37, 35: 37, 72: 37, 94: 37, 101: 37, 5: 36, 92: 36, 26: 35, 64: 35, 79: 35, 84: 35, 85: 35, 116: 35, 18: 34, 24: 34, 82: 34, 44: 33, 81: 33, 100: 33, 32: 32, 86: 32, 90: 32, 25: 31, 43: 31, 80: 31, 112: 31, 47: 29, 57: 29, 104: 29, 38: 28, 50: 28, 54: 28, 66: 28, 69: 27, 98: 27, 119: 27, 123: 27, 132: 27, 70: 26, 113: 26, 115: 26, 62: 25, 96: 25, 17: 24, 56: 24, 59: 24, 105: 24, 117: 24, 88: 23, 129: 23, 19: 22, 33: 22, 95: 22, 108: 22, 53: 21, 67: 21, 68: 21, 114: 21, 118: 21, 133: 21, 10: 20, 51: 20, 125: 19, 40: 18, 42: 18, 37: 17, 103: 17, 14: 16, 106: 15, 130: 15, 65: 14, 107: 13, 109: 13, 61: 12, 127: 9, 135: 9, 134: 8, 124: 7, 120: 5, 52: 3})
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    def tearDown(self):
        reload(dna)

    def test_fit_nonzero_parameters(self):
        c = Counter([1,2,2,3,3,3])
        distribution = model.fit_nonzero_parameters([c],6)
        # Produce output of correct length
        assert_equal(len(distribution), 7)
        # Produce strictly positive parameters
        assert_equal((distribution > 0).all,True)


    def test_log_probability_order(self):
        f = fileinput.input("generated_contigs_test.fna")
        contigs = list(SeqIO.parse(f,"fasta"))
        f.close()
        dna_c1g1 = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c2g1 = dna.DNA(id = c[1].id, seq = str(c[1].seq))
        dna_c3g1 = dna.DNA(id = c[2].id, seq = str(c[2].seq))

        dna_c1g2 = dna.DNA(id = c[-3].id, seq = str(c[-3].seq))
        dna_c2g2 = dna.DNA(id = c[-2].id, seq = str(c[-2].seq))
        dna_c3g2 = dna.DNA(id = c[-1].id, seq = str(c[-1].seq))
        
        cluster1 = [dna_c1g1,dna_c2g1,dna_c3g1]
        cluster2 = [dna_c1g2,dna_c2g2,dna_c3g2]
        for contig in cluster1 + cluster2:
            contig.calculate_signature()

        parameters1 = model.fit_nonzero_parameters([c.signature() for c in cluster1])
        parameters2 = model.fit_nonzero_parameters([c.signature() for c in cluster2])

        s = dna_c1g1.signature()
        log_prob1 = model.log_probability(s,parameters1)
        log_prob2 = model.log_probability(s,parameters2)
        assert_equal(log_prob1>log_prob2,True)
    
    def test_log_probability_basic(self):
        f = fileinput.input("generated_contigs_test.fna")
        contigs = list(SeqIO.parse(f,"fasta"))
        f.close()
        dna_c1g1 = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c2g1 = dna.DNA(id = c[1].id, seq = str(c[1].seq))
        dna_c3g1 = dna.DNA(id = c[2].id, seq = str(c[2].seq))

        dna_c1g2 = dna.DNA(id = c[-3].id, seq = str(c[-3].seq))
        dna_c2g2 = dna.DNA(id = c[-2].id, seq = str(c[-2].seq))
        dna_c3g2 = dna.DNA(id = c[-1].id, seq = str(c[-1].seq))
        
        cluster1 = [dna_c1g1,dna_c2g1,dna_c3g1]
        cluster2 = [dna_c1g2,dna_c2g2,dna_c3g2]
        for contig in cluster1 + cluster2:
            contig.calculate_signature()

        parameters1 = model.fit_nonzero_parameters([c.signature() for c in cluster1])
        parameters2 = model.fit_nonzero_parameters([c.signature() for c in cluster2])

        s1 = dna_c1g1.signature()
        log_prob1 = model.log_probability(s1,parameters1)
        assert_almost_equal(log_prob1, -1000)
        log_prob2 = model.log_probability(s1,parameters2)
        assert_equal(log_prob2,-1000)

        s2 = dna_c1g2.signature()
        log_prob3 = model.log_probability(s2,parameters1)
        assert_equal(log_prob3,-1000)
        log_prob4 = model.log_probability(s2,parameters2)
        assert_equal(log_prob4,-1000)
    
