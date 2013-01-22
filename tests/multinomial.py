#!/usr/bin/env python
from probin.model.composition import multinomial as ml
import random
import fileinput
import unittest
from nose.tools import assert_almost_equal
from Bio import SeqIO

def test_uniform_one_contig_prob():
    f = fileinput.input("data/bambus2.scaffold.linear.fasta.one_contig")
    c = list(SeqIO.parse(f,"fasta"))
    f.close
    signatures = ml.calculate_signatures(4, c)
    k = 4**4
    uniform_prob = [1.0/k]*k
    s=signatures[0]
    log_prob = ml.log_probability(s,uniform_prob)
    assert_almost_equal(log_prob, -1001.20751357)
