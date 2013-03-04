#!/usr/bin/env python
from probin import dna
from collections import Counter
from nose.tools import assert_almost_equal, assert_equal

# testing attribute: signature
class TestDNA(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    def tearDown(self):
        reload(dna)
        
    def test_signature_1(self):
        
        a = dna.DNA(id="ADF",seq="ACTTNACTT")
        correct_signature = Counter([9, 9])
        assert_equal(a.signature, correct_signature)
    
    def test_signature_2(self):
        
        b = dna.DNA(id="ADFA",seq="ACTTTAAACCC")
        true_sign = [9,2, 58, 70, 58, 3, 15, 57]
        correct_signature = Counter(true_sign)
        assert_equal(b.signature, correct_signature)

    def test_signature_3(self):
        
        c = dna.DNA(id="ADADAD",seq='AAAATTTTACGTAGAGCCATTGAGACCTT')
        n_sum = sum(c.signature.values())
        n_len = len("".join(c.seq))
        correct_signature = Counter({21: 2, 85: 2, 228: 2, 5: 1, 136: 1, 144: 1, 157: 1, 167: 1, 173: 1, 183: 1, 57: 1, 188: 1, 61: 1, 193: 1, 197: 1, 211: 1, 84: 1, 220: 1, 221: 1, 112: 1, 241: 1, 116: 1, 245: 1})
        assert_equal(c.signature, correct_signature)
        assert_equal(n_len, n_sum+3)
    
    def test_empty_seq(self):
        
        c = dna.DNA(id="ADADAD",seq='')
        assert_equal(bool(c),True)
         