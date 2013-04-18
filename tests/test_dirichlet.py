#!/usr/bin/env python
from probin.model.composition import dirichlet as model
from probin import dna
import fileinput
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from collections import Counter
import numpy as np
import os
from math import log

cur_dir = os.path.dirname(__file__)

class TestDirichlet(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    def tearDown(self):
        reload(dna)

    def test_fit_nonzero_parameters(self):
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')
        c.calculate_signature()
        d = dna.DNA(id="ADADAD",seq='AAACCCCTATC')
        d.calculate_signature()
        
        alpha_fit = model.fit_nonzero_parameters([c,d])
        # Produce output of correct length
        assert_equal(len(alpha_fit), 136)
        # Produce strictly positive parameters
        assert_equal((alpha_fit > 0).all(),True)

    def test_log_probability_without_fit(self):
        # Example from figure:
        # http://en.wikipedia.org/wiki/File:Beta-binomial_distribution_pmf.png
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')

        signature = np.array([6,4])

        alpha = np.array([600,400])
        p = model.log_probability_test(signature,alpha)
        # 10 choose 6 = 210
        assert_almost_equal(p + np.log(210),-1.38799, places=3)


        signature_m = np.zeros((1,2))
        signature_m[0,:] = np.array([6,4])

        p2 = model.neg_log_probability_l(alpha,signature_m,1,10)
        assert_almost_equal(-p2+np.log(210),-1.38799, places=3)


    def test_log_probability_list(self):
        c = dna.DNA(id="ADADAD", seq='ACTTTAAACCC')
        c.calculate_signature()
        d = dna.DNA(id="ADADAD", seq='ACTTTACGAACCC')
        d.calculate_signature()
        dna_l = [c,d]
        alpha = [3.0]*c.kmer_hash_count

        sig_mat = np.zeros((2,c.kmer_hash_count))
        for key,cnt in c.signature.iteritems():
            sig_mat[0,key] += cnt

        for key,cnt in d.signature.iteritems():
            sig_mat[1,key] += cnt

        p = model.neg_log_probability_l(alpha,sig_mat,2,np.array([10,8]))
        assert_almost_equal(p,88.31776,places=4)
        

    def test_log_probability_order(self):
        file_name = os.path.join(cur_dir,"..","data/generated_contigs_test.fna")
        f = fileinput.input(file_name)
        c = list(SeqIO.parse(f,"fasta"))
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

        parameters1 = model.fit_nonzero_parameters(cluster1)
        parameters2 = model.fit_nonzero_parameters(cluster2)

        log_prob1 = model.log_probability(dna_c1g1,parameters1)
        log_prob2 = model.log_probability(dna_c1g1,parameters2)
        assert_equal(log_prob1>log_prob2,True)
    
    def test_log_probability_full(self):
        file_name = os.path.join(cur_dir,"..","data/generated_contigs_test.fna")
        f = fileinput.input(file_name)
        c = list(SeqIO.parse(f,"fasta"))
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

        parameters1 = model.fit_nonzero_parameters(cluster1)
        parameters2 = model.fit_nonzero_parameters(cluster2)

        # These testa are probably too shaky, due to the
        # numerical optimization for finding the parameters
        s1 = dna_c1g1
        log_prob1 = model.log_probability(s1,parameters1)
        assert_almost_equal(log_prob1/10000.0, -0.450, places = 1)
        log_prob2 = model.log_probability(s1,parameters2)
        assert_almost_equal(log_prob2/10000.0,-0.4676, places = 2)

        s2 = dna_c1g2
        log_prob3 = model.log_probability(s2,parameters1)
        assert_almost_equal(log_prob3/10000.0,-0.517, places = 2)
        log_prob4 = model.log_probability(s2,parameters2)
        assert_almost_equal(log_prob4/10000.0,-0.483, places = 2)
