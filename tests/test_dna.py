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
    
    def test_pseudo_count1(self):
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')
        c.calculate_signature()
        n_sum = sum(c.signature.values())
        n_len = len("".join(c.seq))
        true_sign = [9,2, 58, 70, 58, 3, 15, 57]
        assert_equal(c.pseudo_count(9),2 )
        assert_equal(c.pseudo_count(2),2 )
        assert_equal(c.pseudo_count(58),3 )
        assert_equal(c.pseudo_count(70),2 )
        assert_equal(c.pseudo_count(3),2 )
        assert_equal(c.pseudo_count(15),2 )
        assert_equal(c.pseudo_count(57),2 )
        
        assert_equal(c.pseudo_count(10), 1)

    def test_pseudo_counts(self):
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')
        c.calculate_signature()
        true_sign = [9,2, 58, 70, 58, 3, 15, 57]
        pseudo_sign = [i for i in c.pseudo_counts]
        assert_equal(len(pseudo_sign), 136)
        assert_equal(pseudo_sign[9],2)
        
    def test_pseudo_count2(self):
        a = dna.DNA(id="ADF",seq="ACTTNACTT")
        a.calculate_signature()
        correct_signature = Counter([9, 9])
        assert_equal(a.pseudo_count(9) == 3, True)
        assert_equal(a.pseudo_count(10) == 1,True)
        

    def test_empty_seq(self):
        c = dna.DNA(id="ADADAD",seq='')
        assert_equal(bool(c),True)
         
    
    def test_signature_calculation_is_not_in_constructor(self):
        a = dna.DNA(id="ADADAD",seq='AAAATTTTACGTAGAGCCATTGAGACCTT')
        assert_is_none(a.signature)
        
        
    def test_split_seq_to_signatures1(self):
        test_seq = 'AAAATTTTACGTAGAGCCATTGAGACCTT'
        a = dna.DNA(id="ADADAD",seq=test_seq)
        sigs = a.split_seq_to_signatures(10,3)
        a.calculate_signature()
        original_sig = a.signature

        assert_equal(len(sigs), 3)
        # All kmers should be in the test_seq.
        assert_equal(all([key in original_sig for key in sigs[0].keys()]), True)
        assert_equal(all([key in original_sig for key in sigs[1].keys()]), True)       
        assert_equal(all([key in original_sig for key in sigs[2].keys()]), True)

