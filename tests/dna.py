#!/usr/bin/env python
from probin.dna import DNA
from collections import Counter
import random
import fileinput
import unittest
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from numpy import array

# testing attribute: signature
def test_signature_1():
    DNA.generate_kmer_hash(4)
    a = DNA(id="ADF",seq="ACTTNACTT")
    correct_signature = Counter([53, 53])
    assert_equal(a.signature, correct_signature)

def test_signature_2():
    b = DNA(id="ADFA",seq="ACTTTAAACCCACACACAACATTTGGAAAGGAGAGAGCCATTA")
    correct_signature = Counter({204: 3, 136: 2, 153: 2, 197: 2, 84: 2, 213: 2, 221: 2, 128: 1, 149: 1, 160: 1, 165: 1, 169: 1, 173: 1, 49: 1, 53: 1, 183: 1, 188: 1, 192: 1, 65: 1, 195: 1, 80: 1, 89: 1, 223: 1, 101: 1, 252: 1, 240: 1, 241: 1, 243: 1, 245: 1, 124: 1, 125: 1, 21: 1})
    assert_equal(b.signature, correct_signature)

def test_signature_3():
    c = DNA(id="ADADAD",seq='AAAATTTTACGTAGAGCCATTGAGACCTT')
    n_sum = sum(c.signature.values())
    n_len = len("".join(c.seq))
    correct_signature = Counter({21: 2, 85: 2, 228: 2, 5: 1, 136: 1, 144: 1, 157: 1, 167: 1, 173: 1, 183: 1, 57: 1, 188: 1, 61: 1, 193: 1, 197: 1, 211: 1, 84: 1, 220: 1, 221: 1, 112: 1, 241: 1, 116: 1, 245: 1})
    assert_equal(c.signature, correct_signature)
    assert_equal(n_len, n_sum+3)
def test_empty_seq():
    c = DNA(id="ADADAD",seq='')
    assert_equal(bool(c),True)
    
