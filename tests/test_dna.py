#!/usr/bin/env python
from probin import dna
from collections import Counter
from nose.tools import assert_almost_equal, assert_equal, assert_is_none

# testing attribute: signature
class TestDNA(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    def tearDown(self):
        reload(dna)
        
    def test_signature_1(self):
        a = dna.DNA(id="ADF",seq="ACTTNACTT")
        a.calculate_signature()
        correct_signature = Counter([9, 9])
        assert_equal(a.signature, correct_signature)
    
    def test_signature_2(self):
        b = dna.DNA(id="ADFA",seq="ACTTTAAACCC")
        b.calculate_signature()
        true_sign = [9,2, 58, 70, 58, 3, 15, 57]
        correct_signature = Counter(true_sign)
        assert_equal(b.signature, correct_signature)

    def test_signature_3(self):
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')
        c.calculate_signature()
        n_sum = sum(c.signature.values())
        n_len = len("".join(c.seq))
        true_sign = [9,2, 58, 70, 58, 3, 15, 57]
        correct_signature = Counter(true_sign)
        assert_equal(c.signature, correct_signature)
        assert_equal(n_len, n_sum+3)
    
    def test_empty_seq(self):
        c = dna.DNA(id="ADADAD",seq='')
        assert_equal(bool(c),True)
         
    
    def test_signature_calculation_is_not_in_constructor(self):
        a = dna.DNA(id="ADADAD",seq='AAAATTTTACGTAGAGCCATTGAGACCTT')
        assert_is_none(a.signature)
