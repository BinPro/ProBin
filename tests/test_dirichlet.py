#!/usr/bin/env python
from probin.model.composition import dirichlet as model
from probin import dna
import fileinput
from nose.tools import assert_almost_equal, assert_equal
from Bio import SeqIO
from collections import Counter
import numpy as np
from math import log

class TestDirichlet(object):
    def setUp(self):
        dna.DNA.generate_kmer_hash(4)
    def tearDown(self):
        reload(dna)

    def test_fit_nonzero_parameters(self):
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')
        c.calculate_signature()
        
        alpha_fit = model.fit_nonzero_parameters([c],4)
        # Produce output of correct length
        assert_equal(len(alpha_fit), 4)
        # Produce strictly positive parameters
        assert_equal((alpha_fit > 0).all(),True)

    def test_log_probability_without_fit(self):
        c = dna.DNA(id="ADADAD",seq='ACTTTAAACCC')
        c.calculate_signature()

        alpha = np.array([600,400])
        p = model.log_probability(c,alpha)
        assert_almost_equal(p,-992.1644316)


    def test_log_probability_order(self):
        f = fileinput.input("generated_contigs_test.fna")
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

        parameters1 = model.fit_nonzero_parameters(cluster1,dna.DNA.kmer_hash_count)
        parameters2 = model.fit_nonzero_parameters(cluster2,dna.DNA.kmer_hash_count)

        log_prob1 = model.log_probability(dna_c1g1,parameters1)
        log_prob2 = model.log_probability(dna_c1g1,parameters2)
        assert_equal(log_prob1>log_prob2,True)
    
    def test_log_probability_full(self):
        f = fileinput.input("generated_contigs_test.fna")
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

        parameters1 = model.fit_nonzero_parameters(cluster1,dna.DNA.kmer_hash_count)
        parameters2 = model.fit_nonzero_parameters(cluster2,dna.DNA.kmer_hash_count)

        # These testa are probably too shaky, due to the
        # numerical optimization for finding the parameters
        s1 = dna_c1g1
        log_prob1 = model.log_probability(s1,parameters1)
        assert_almost_equal(log_prob1/10000.0, -0.52, places = 2)
        log_prob2 = model.log_probability(s1,parameters2)
        assert_almost_equal(log_prob2/10000.0,-0.5358, places = 2)

        s2 = dna_c1g2
        log_prob3 = model.log_probability(s2,parameters1)
        assert_almost_equal(log_prob3/10000.0,-0.587, places = 2)
        log_prob4 = model.log_probability(s2,parameters2)
        assert_almost_equal(log_prob4/10000.0,-0.5591, places = 2)
    
