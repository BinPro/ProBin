#!/usr/bin/env python
from probin.model.composition import multinomial as ml
from probin.dna import DNA
import random
import fileinput
import unittest
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from numpy import array
from collections import Counter

# testing function: log_probability
def test_uniform_one_contig_prob():
    DNA.generate_kmer_hash(4)
    f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
    c = list(SeqIO.parse(f,"fasta"))
    f.close()
    dna_c = DNA(id = c[0].id, seq = str(c[0].seq))
    s = dna_c.signature
    k = 4**4
    uniform_prob = {}
    for i,cnt in s.items():
        uniform_prob[i] = 1./k
    log_prob = ml.log_probability(s,uniform_prob)
    print log_prob
    assert_almost_equal(log_prob, -3791.05738056)

# testing function: calculate_signatures
def test_signatures_one_contig_basic():
    f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
    c = list(SeqIO.parse(f,"fasta"))
    f.close()
    dna_c = DNA(id = c[0].id, seq = str(c[0].seq))
    calculated_signature = dna_c.signature
    correct_signature = CORRECT_SIGNATURES_ONE_CONTIG
    assert_equal(calculated_signature,correct_signature)

# testing function: fit_parameters
def test_parameters_one_contig_basic():
    f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
    c = list(SeqIO.parse(f,"fasta"))
    f.close()
    dna_c = DNA(id = c[0].id, seq = str(c[0].seq))
    c_sig = CORRECT_SIGNATURES_ONE_CONTIG
    n = sum(c_sig.values())
    correct_parameters = {}
    for i,v in c_sig.items():
        correct_parameters[i] = v/float(n)
    calculated_parameters = ml.fit_parameters(dna_c.signature)
    assert_equal(calculated_parameters, correct_parameters)

def test_signaturs_large_genome():
    f = fileinput.input("data/8M_genome.fna")
    c= list(SeqIO.parse(f,"fasta"))
    f.close
    dna_c = DNA(id = c[0].id, seq = str(c[0].seq))
    calculated_parameters = ml.fit_parameters(dna_c.signature)
    assert_equal(len(calculated_parameters), 136)



## "Constants"
CORRECT_SIGNATURES_ONE_CONTIG = Counter({21: 73, 85: 73, 97: 64, 192: 62, 112: 60, 84: 59, 165: 59, 193: 59, 133: 58, 101: 56, 128: 56, 176: 56, 199: 56, 149: 55, 37: 53, 105: 53, 89: 52, 132: 52, 69: 51, 93: 51, 65: 49, 73: 49, 81: 49, 181: 49, 213: 49, 33: 48, 96: 48, 109: 48, 203: 47, 204: 47, 214: 47, 148: 46, 195: 46, 201: 46, 216: 46, 240: 46, 243: 46, 53: 45, 225: 45, 184: 44, 188: 44, 217: 44, 161: 43, 219: 43, 168: 42, 200: 42, 242: 41, 141: 40, 208: 40, 117: 39, 236: 39, 246: 39, 253: 39, 77: 38, 129: 38, 197: 38, 113: 37, 136: 37, 144: 37, 151: 37, 205: 37, 5: 36, 152: 36, 116: 35, 124: 35, 140: 35, 177: 35, 207: 35, 224: 35, 49: 34, 156: 34, 209: 34, 173: 33, 215: 33, 220: 33, 172: 32, 221: 32, 252: 32, 160: 31, 179: 31, 237: 31, 241: 31, 153: 29, 169: 29, 247: 29, 61: 28, 137: 28, 180: 28, 185: 28, 164: 27, 226: 27, 232: 27, 239: 27, 251: 27, 80: 26, 223: 26, 227: 26, 196: 25, 248: 25, 17: 24, 183: 24, 212: 24, 233: 24, 255: 24, 108: 23, 230: 23, 120: 22, 145: 22, 157: 22, 211: 22, 100: 21, 163: 21, 191: 21, 228: 21, 249: 21, 254: 21, 121: 20, 245: 20, 235: 19, 45: 18, 189: 18, 125: 17, 135: 17, 229: 16, 210: 15, 231: 15, 244: 14, 147: 13, 167: 13, 68: 12, 198: 9, 250: 9, 238: 8, 187: 7, 175: 5, 57: 3})
