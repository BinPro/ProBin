#!/usr/bin/env python
from Bio import SeqIO
from random import randint 

def sample_contigs(genome, n, min_length, max_length):
    contigs = []
    for i in range(n):
        l = randint(min_length, max_length)
        start = randint(0, (len(genome)-l))
        contigs.append(str(genome.seq[start:start+l]))
    return contigs
