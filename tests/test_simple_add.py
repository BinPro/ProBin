#!/usr/bin/env python
import os
import sys
import fileinput
from nose.tools import assert_equal, assert_true, assert_almost_equal
import numpy as np
import scipy.stats
from collections import Counter
from Bio import SeqIO

from probin import dna
from probin.model.coverage import isotropic_gaussian as ig
from probin.model.composition import multinomial as ml
from probin.model.combined import simple_add as model

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))

class TestSimpleAdd(object):
    ## "Constants"
    CORRECT_SIGNATURES_ONE_CONTIG = Counter({0: 73, 1: 73, 27: 64, 73: 62, 75: 60, 15: 59, 21: 59, 58: 59, 7: 58, 12: 56, 74: 56, 77: 56, 102: 56, 3: 55, 13: 53, 55: 53, 45: 52, 63: 52, 4: 51, 31: 51, 2: 49, 11: 49, 16: 49, 20: 49, 48: 49, 28: 48, 41: 48, 78: 48, 83: 47, 122: 47, 126: 47, 49: 46, 60: 46, 76: 46, 91: 46, 110: 46, 111: 46, 9: 45, 29: 45, 46: 44, 87: 44, 97: 44, 30: 43, 121: 43, 93: 42, 99: 42, 131: 41, 36: 40, 71: 40, 8: 39, 39: 39, 89: 39, 128: 39, 6: 38, 22: 38, 34: 38, 23: 37, 35: 37, 72: 37, 94: 37, 101: 37, 5: 36, 92: 36, 26: 35, 64: 35, 79: 35, 84: 35, 85: 35, 116: 35, 18: 34, 24: 34, 82: 34, 44: 33, 81: 33, 100: 33, 32: 32, 86: 32, 90: 32, 25: 31, 43: 31, 80: 31, 112: 31, 47: 29, 57: 29, 104: 29, 38: 28, 50: 28, 54: 28, 66: 28, 69: 27, 98: 27, 119: 27, 123: 27, 132: 27, 70: 26, 113: 26, 115: 26, 62: 25, 96: 25, 17: 24, 56: 24, 59: 24, 105: 24, 117: 24, 88: 23, 129: 23, 19: 22, 33: 22, 95: 22, 108: 22, 53: 21, 67: 21, 68: 21, 114: 21, 118: 21, 133: 21, 10: 20, 51: 20, 125: 19, 40: 18, 42: 18, 37: 17, 103: 17, 14: 16, 106: 15, 130: 15, 65: 14, 107: 13, 109: 13, 61: 12, 127: 9, 135: 9, 134: 8, 124: 7, 120: 5, 52: 3})
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)

    def tearDown(self):
        reload(dna)
            
    def test_log_probability(self):
        f = fileinput.input(os.path.join(data_path,"bambus2.scaffold.linear.fasta.one_contig"))
        c = list(SeqIO.parse(f,"fasta"))
        f.close()
        
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()
        k = 4**4
        uniform_prob = np.ones((dna.DNA.kmer_hash_count))
        for i,cnt in dna_c.signature.items():
            uniform_prob[i] = 1./k
        log_prob = ml.log_probability(dna_c,uniform_prob)

        mu = np.log(np.array([0.5,3.0,5.0]))
        sigma = 0.5
        x = np.log(np.array([0.5,3.0,5.0]))
        cov_matrix = np.array([x])

        p_ig = ig.log_pdf(x,mu,sigma)
        p_test = model.log_probability(dna_c,cov_matrix,uniform_prob,mu,sigma)
        assert_equal(p_ig+log_prob,p_test)

        
    def test_fit_parameters_single_cluster(self):
        # A sample with two contigs, each with three data points.
        
        f = fileinput.input(os.path.join(data_path,"bambus2.scaffold.linear.fasta.one_contig"))
        c = list(SeqIO.parse(f,"fasta"))
        f.close()
        
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()

        x = np.log(np.array([[0.5,3.0,2.0],[3.0,1.0,1.0]]))
        cov_matrix = x
        # one cluster, two contigs
        exp_clust = np.array([[1.0],[1.0]])
        mu0 = np.log(np.array([0.5,3.0])).sum()/2.0
        mu1 = np.log([3.0,1.0]).sum()/2.0
        mu2 = np.log([2.0,1.0]).sum()/2.0
        correct_mu = np.array([mu0,mu1,mu2])

        diff_vec = [(x[i,0] - mu0)**2 + (x[i,1] -mu1)**2 + (x[i,2]-mu2)**2 for i in range(2)]
        correct_sigma = np.array(diff_vec).sum()

        correct_sigma /= 2.0
        
        c_sig = self.CORRECT_SIGNATURES_ONE_CONTIG
        n = sum(c_sig.values())
        correct_parameters_mul = np.zeros(dna.DNA.kmer_hash_count)
        for i,v in c_sig.iteritems():
            correct_parameters_mul[i] += v*2 + 1
        correct_parameters_mul/=np.sum(correct_parameters_mul)
        cal_prob_v,cal_mu,cal_sigma = model.fit_nonzero_parameters([dna_c,dna_c],cov_matrix,expected_clustering=exp_clust)

        assert_equal((cal_prob_v==correct_parameters_mul).all(),True)
        assert_equal((cal_mu == correct_mu).all(),True)
        assert_equal((cal_sigma == correct_sigma).all(),True)



    def test_fit_parameters_two_clusters(self):
        # A sample with four contigs, each with three data points.

        f = fileinput.input(os.path.join(data_path,"bambus2.scaffold.linear.fasta.one_contig"))
        c = list(SeqIO.parse(f,"fasta"))
        f.close()
        
        dna_c = dna.DNA(id = c[0].id, seq = str(c[0].seq))
        dna_c.calculate_signature()
        
        x = np.log(np.array([[0.5,3.0,2.0],[3.0,1.0,1.0],
                             [2.5,1.5,1.0],[1.0,1.0,2.0]]))
        cov_matrix = x
        # two clusters, four contigs
        exp_clust = np.array([[0.7,0.3],
                              [0.1,0.9],
                              [0.5,0.5],
                              [0.2,0.8]])

        mu00 = (np.log(0.5)*0.7+np.log(3.0)*0.1+np.log(2.5)*0.5+np.log(1.0)*0.2)/1.5
        mu10 = (np.log(3.0)*0.7+np.log(1.0)*0.1+np.log(1.5)*0.5+np.log(1.0)*0.2)/1.5
        mu20 = (np.log(2.0)*0.7+np.log(1.0)*0.1+np.log(1.0)*0.5+np.log(2.0)*0.2)/1.5

        mu01 = (np.log(0.5)*0.3+np.log(3.0)*0.9+np.log(2.5)*0.5+np.log(1.0)*0.8)/2.5
        mu11 = (np.log(3.0)*0.3+np.log(1.0)*0.9+np.log(1.5)*0.5+np.log(1.0)*0.8)/2.5
        mu21 = (np.log(2.0)*0.3+np.log(1.0)*0.9+np.log(1.0)*0.5+np.log(2.0)*0.8)/2.5

        correct_mu = np.array([[mu00,mu10,mu20],[mu01,mu11,mu21]])
        mu = correct_mu

        sigma_test0 = np.array([((x[i,0] - mu[0,0])**2 + (x[i,1] -mu[0,1])**2 + (x[i,2]-mu[0,2])**2)*exp_clust[i,0] for i in range(4)]).sum()
        sigma_test0 /= 1.5

        sigma_test1 = np.array([((x[i,0] - mu[1,0])**2 + (x[i,1] -mu[1,1])**2 + (x[i,2]-mu[1,2])**2)*exp_clust[i,1] for i in range(4)]).sum()
        sigma_test1 /= 2.5
        
        correct_sigma = np.array([sigma_test0,sigma_test1])


        c_sig = self.CORRECT_SIGNATURES_ONE_CONTIG
        n = sum(c_sig.values())
        correct_parameters_mul0 = np.zeros(dna.DNA.kmer_hash_count)
        correct_parameters_mul1 = np.zeros(dna.DNA.kmer_hash_count)
        for i,v in c_sig.iteritems():
            correct_parameters_mul0[i] += v*(0.7+0.1+0.5+0.2) + 1
            correct_parameters_mul1[i] += v*(0.3+0.9+0.5+0.8) + 1
        correct_parameters_mul0/=np.sum(correct_parameters_mul0)
        correct_parameters_mul1/=np.sum(correct_parameters_mul1)
        correct_parameters_mul = np.array([correct_parameters_mul0,
                                           correct_parameters_mul1])
        cal_prob_v,cal_mu,cal_sigma = model.fit_nonzero_parameters([dna_c,dna_c,dna_c,dna_c],cov_matrix,expected_clustering=exp_clust)

        assert_equal((np.abs(cal_prob_v-correct_parameters_mul)<1e-5).all(),True)
                
        assert_equal((np.abs(cal_mu - correct_mu)<1e-7).all(),True)
        assert_equal((np.abs(cal_sigma - correct_sigma)<1e-7).all(),True)
