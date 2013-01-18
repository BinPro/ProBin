#!/usr/bin/env python
from probin.helpers.sequence_signature import SequenceSignature as ss
from  itertools import product
from numpy import array, log
from scipy.misc import factorial
from collections import defaultdict

def calculate_signatures(kmer_length,contigs):
    signatures = []
    kmers = possible_kmers(kmer_length)
    for c in contigs:
        signature = ss(kmer_length,c, kmers).kmer_frequencies
        signatures.append(signature[1:])
    return signatures

def possible_kmers(k):
    kmer_dict = defaultdict(int)
    kmer_list = [''.join(x) for x in product('ATGC', repeat=k)]
    kmer_list.insert(0,'REST')
    for i in range(len(kmer_list)):
        kmer_dict[kmer_list[i]] = i
    return kmer_dict

def probability(signature, prob_vector):
    phi = sum(signature)
    prod = factorial(phi)
    for j in range(len(prob_vector)):
        denom = factorial(signature[j])
        prod *= (prob_vector[j]**signature[j])/denom
    return prod

def log_probability(signature, prob_vector):
    phi = sum(signature)
    log_prod = 0
    for j in range(len(prob_vector)):
        denom = log_fac(signature[j])
        log_prod += (log(prob_vector[j])*signature[j]) - denom
    return log_prod + log_fac(phi)

def log_fac(i):
    return sum(log(range(1,i+1)))
