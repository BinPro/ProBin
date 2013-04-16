#!/usr/bin/env python
from probin import dna
from collections import Counter
from nose.tools import assert_almost_equal, assert_equal, assert_is_none, raises

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

    def test_pseudo_count2(self):
        a = dna.DNA(id="ADF",seq="ACTTNACTT")
        a.calculate_signature()
        correct_signature = Counter([9, 9])
        assert_equal(a.pseudo_count(9) == 3, True)
        assert_equal(a.pseudo_count(10) == 1,True)

    def test_pseudo_counts(self):
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')
        c.calculate_signature()
        true_sign = [9,2, 58, 70, 58, 3, 15, 57]
        pseudo_sign = [i for i in c.pseudo_counts]
        assert_equal(len(pseudo_sign), 136)
        assert_equal(pseudo_sign[9],2)
        assert_equal(sum(c.pseudo_counts),len(true_sign)+dna.DNA.kmer_hash_count)
        

    def test_empty_seq(self):
        c = dna.DNA(id="ADADAD",seq='')
        assert_equal(c.full_seq,'')
        assert_equal(len(c),0)
        assert_equal(bool(c),False)
    
    @raises(Exception)
    def test_changing_seq(self):
        c = dna.DNA(id="ADADAD",seq='')
        c.seq = "ACTTTAAACCC"
         
    
    def test_signature_calculation_is_not_in_constructor(self):
        a = dna.DNA(id="ADADAD",seq='AAAATTTTACGTAGAGCCATTGAGACCTT')
        assert_is_none(a.signature)

    def test_pick_part(self):
        test_seq = 'AAAATTTTACGTAGAGCCATTGAGACCTT'
        a = dna.DNA(id="ADADAD",seq=test_seq)
        part, left_right = a._pick_part(a.full_seq,0,len(test_seq),5)
        left = left_right[0]
        right = left_right[1]
        assert_equal(len(part[0]),5)
        assert_equal(left[0] + part[0] + right[0], test_seq)
        assert_equal(len(left[0]) + len(right[0]), len(test_seq)-5)
        # Start position
        assert_equal(left[1],0)
        assert_equal(test_seq[part[1]:part[1]+5],part[0])
        assert_equal(len(test_seq[right[1]:]),len(right[0]))
        assert_equal(test_seq[right[1]:],right[0])


    def test_split_seq_random(self):
        test_seq = 'AAAATTTTACGTAGAGCCATTGAGACCTT'
        a = dna.DNA(id="ADADAD",seq=test_seq)
        seqs = a.split_seq_random(10,3)

        assert_equal(len(seqs), 3)
        assert_equal(len(seqs[0].full_seq),10)

        # All sequences should be in the test_seq.
        assert_equal((seqs[0].full_seq in a.full_seq), True)
        assert_equal((seqs[1].full_seq in a.full_seq), True)
        assert_equal((seqs[2].full_seq in a.full_seq), True)

        assert_equal(seqs[0].start_position>=0,True)
        assert_equal(test_seq[seqs[0].start_position:][:10],seqs[0].full_seq)

    def test_split_seq(self):
        test_seq = 'AAAATTTTACGTAGAGCCATTGAGACCTT'
        a = dna.DNA(id="ADADAD",seq=test_seq)
        seqs = a.split_seq(10)

        assert_equal(len(seqs), 3)
        assert_equal(len(seqs[0].full_seq),10)
        assert_equal(len(seqs[1].full_seq),10)
        assert_equal(len(seqs[2].full_seq),9)

        assert_equal(seqs[0].full_seq+seqs[1].full_seq+seqs[2].full_seq, test_seq)

        assert_equal(seqs[0].start_position==0,True)
        assert_equal(test_seq[seqs[0].start_position:][:10],seqs[0].full_seq)

    def test_len_seq(self):
        test_seq = 'AAAATTTTACGTAGAGCCATTGAGACCTT'
        a = dna.DNA(id='ADADAD',seq=test_seq)
        assert_equal(len(a),len(test_seq))
        assert_equal(a._sequence_length,len(a))
