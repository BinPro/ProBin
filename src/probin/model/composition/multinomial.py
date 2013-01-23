#!/usr/bin/env python
from probin.helpers.sequence_signature import SequenceSignature as SS
from probin.helpers.sequence_signature import possible_kmers
from probin.helpers.misc import log_fac
from numpy import log


def calculate_signatures(kmer_length,contigs):
    signatures = []
    kmers = possible_kmers(kmer_length)
    for c in contigs:
        signature = SS(kmer_length,c, kmers).kmer_frequencies
        signatures.append(signature[1:])
    return signatures

def fit_parameters(kmer_length, contigs):
    signatures = calculate_signatures(kmer_length, contigs)
    return [[n/float(sum(sig)) for n in sig] for sig in signatures]

def log_probability(signature, prob_vector):
    phi = sum(signature)
    log_prod = 0
    for j in range(len(prob_vector)):
        denom = log_fac(signature[j])
        log_prod += (log(prob_vector[j])*signature[j]) - denom
    return log_prod + log_fac(phi)
