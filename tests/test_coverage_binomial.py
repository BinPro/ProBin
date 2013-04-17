#!/usr/bin/env python
from probin.model.coverage import binomial
from probin import dna
import fileinput
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from collections import Counter
import numpy as np
from numpy import log
import sys

class TestCoverageBionomial(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    def tearDown(self):
        reload(dna)
    # testing function: log_probability
    def test_log_probability(self):
        """ With a q_prime = 0.5 the binomial probability
        should be easy to estimate. 
        """
        c = dna.DNA(id = "ADADAD", seq="A")
        c.mapping_reads = np.array([1])
        lp = binomial.log_probability(c, np.array([0.5]),np.array([2]),1)
        # The probability of getting tail from a single coin toss
        assert_almost_equal(lp,log(0.5))

        # Different L values
        lp = binomial.log_probability(c, np.array([0.5]),np.array([2]),2)
        assert_almost_equal(lp,log(0.375))
        
        # Different sequence length
        c = dna.DNA(id = "ADADAD", seq="AT")
        c.mapping_reads = np.array([1])
        lp = binomial.log_probability(c, np.array([0.5]),np.array([2]),2)
        assert_almost_equal(lp,log(0.5))
        
    def test_log_probability_series(self):
        """ Basic series test """
        c = dna.DNA(id = "ADADAD", seq="A")
        c.mapping_reads = np.array([1,0,2])
        lp = binomial.log_probability(c, np.array([0.5,0.5,0.5]),np.array([2,2,2]),1)

        assert_almost_equal(lp, log(0.5)+log(0.25)+log(0.25))
        
        # Different q-values
        lp  = binomial.log_probability(c, np.array([0.5,0.1,0.9]),np.array([2,2,2]),1)
        assert_almost_equal(lp,log(0.5)+log(0.81)+log(0.81))
        
        # Different M-values
        lp =  binomial.log_probability(c, np.array([0.5,0.5,0.5]),np.array([2,3,4]),1)
        assert_almost_equal(lp,log(0.5)+log(0.125)+log(0.375))


    # Testing fit_nonzero_parameters
    def test_fit_nonzero_parameters(self):
        c = dna.DNA(id = "ADADAD", seq="A")
        c.mapping_reads = np.array([2,2,2])
        # Same M
        M = np.array([4,4,4])
        pseudo_par = binomial.fit_nonzero_parameters([c],M)
        assert_equal((pseudo_par==3/float(5)).all(),True)
        
        # Different M:s
        M = np.array([2,6,10])
        pseudo_par = binomial.fit_nonzero_parameters([c],M)
        true_par = np.array([3/float(3),3/float(7),3/float(11)])
        assert_equal(pseudo_par[0],true_par[0])
        assert_equal(pseudo_par[1],true_par[1])
        assert_equal(pseudo_par[2],true_par[2])
