#!/usr/bin/env python
from Bio import SeqIO
from random import randint 

def sample_contigs(genome, n, min_length, max_length):
    contigs = []
    identifier = ">" + genome.id
    for i in range(n):
        l = randint(min_length, max_length)
        start = randint(0, (len(genome)-l))
        contigs.append(identifier + " contig_number " + str(i) + "\n" + str(genome.seq[start:start+l]))
    return contigs
