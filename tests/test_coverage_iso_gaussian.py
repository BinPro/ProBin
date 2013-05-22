#!/usr/bin/env python
import os
import sys
import fileinput
from nose.tools import assert_equal, assert_true, assert_almost_equal
import numpy as np
import scipy.stats
from Bio import SeqIO

from probin import dna
from probin.binning import kmeans
from probin.binning import statistics
from probin.model.coverage import isotropic_gaussian as model

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))

class TestIsotropicGaussian(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)

    def tearDown(self):
        reload(dna)
            
    def test_pdf(self):
        mu = np.log(np.array([0.5,3.0,5.0]))
        sigma = 0.5
        x = np.log(np.array([0.5,3.0,5.0]))
        p0 = scipy.stats.norm.pdf(x[0],loc=mu[0],scale=sigma)
        p1 = scipy.stats.norm.pdf(x[1],loc=mu[1],scale=sigma)
        p2 = scipy.stats.norm.pdf(x[2],loc=mu[2],scale=sigma)
        p_test = model.pdf(x,mu,sigma)
        assert_equal(p0*p1*p2,p_test)
    def test_log_pdf(self):
        mu = np.log(np.array([0.5,3.0,5.0]))
        sigma = 0.5
        x = np.log(np.array([0.5,3.0,5.0]))
        p = model.pdf(x,mu,sigma)
        p_test = model.log_pdf(x,mu,sigma)

        assert_almost_equal(np.log(p),p_test)

    def test_fit_parameters_single_cluster(self):
        # A sample with two contigs, each with three data points.
        
        x = np.log(np.array([[0.5,3.0,2.0],[3.0,1.0,1.0]]))
        # one cluster, two contigs
        exp_clust = np.array([[1.0],[1.0]])
        mu0 = np.log(np.array([0.5,3.0])).sum()/2.0
        mu1 = np.log([3.0,1.0]).sum()/2.0
        mu2 = np.log([2.0,1.0]).sum()/2.0
        mu,sigma = model.fit_parameters(x,expected_clustering = exp_clust)
        assert_equal(mu[0,0],mu0)
        assert_equal(mu[1,0],mu1)
        assert_equal(mu[2,0],mu2)

        sigma_test = np.array([(x[i,0] - mu[0])**2 + (x[i,1] -mu[1])**2 + (x[i,2]-mu[2])**2 for i in range(2)]).sum()
        sigma_test /= 2.0
        
        assert_equal(sigma[0],sigma_test)


    def test_fit_parameters_two_clusters(self):
        # A sample with four contigs, each with three data points.
        
        x = np.log(np.array([[0.5,3.0,2.0],[3.0,1.0,1.0],
                             [2.5,1.5,1.0],[1.0,1.0,2.0]]))
        # two clusters, four contigs
        exp_clust = np.array([[0.7,0.3],[0.1,0.9],
                              [0.5,0.5],[0.2,0.8]])

        mu00 = (np.log(0.5)*0.7+np.log(3.0)*0.1+np.log(2.5)*0.5+np.log(1.0)*0.2)/1.5
        mu10 = (np.log(3.0)*0.7+np.log(1.0)*0.1+np.log(1.5)*0.5+np.log(1.0)*0.2)/1.5
        mu20 = (np.log(2.0)*0.7+np.log(1.0)*0.1+np.log(1.0)*0.5+np.log(2.0)*0.2)/1.5

        mu01 = (np.log(0.5)*0.3+np.log(3.0)*0.9+np.log(2.5)*0.5+np.log(1.0)*0.8)/2.5
        mu11 = (np.log(3.0)*0.3+np.log(1.0)*0.9+np.log(1.5)*0.5+np.log(1.0)*0.8)/2.5
        mu21 = (np.log(2.0)*0.3+np.log(1.0)*0.9+np.log(1.0)*0.5+np.log(2.0)*0.8)/2.5
        mu,sigma = model.fit_parameters(x,expected_clustering = exp_clust)
        assert_almost_equal(mu[0,0],mu00)
        assert_almost_equal(mu[1,0],mu10)
        assert_almost_equal(mu[2,0],mu20)

        assert_almost_equal(mu[0,1],mu01)
        assert_almost_equal(mu[1,1],mu11)
        assert_almost_equal(mu[2,1],mu21)

        sigma_test0 = np.array([((x[i,0] - mu[0,0])**2 + (x[i,1] -mu[1,0])**2 + (x[i,2]-mu[2,0])**2)*exp_clust[i,0] for i in range(4)]).sum()
        sigma_test0 /= 1.5

        sigma_test1 = np.array([((x[i,0] - mu[0,1])**2 + (x[i,1] -mu[1,1])**2 + (x[i,2]-mu[2,1])**2)*exp_clust[i,1] for i in range(4)]).sum()
        sigma_test1 /= 2.5
        
        assert_almost_equal(sigma[0],sigma_test0)
        assert_almost_equal(sigma[1],sigma_test1)
