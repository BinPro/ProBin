#!/usr/bin/env python
from ProBin.Helpers.sequence_signature import Sequence_signature as ss
from  itertools import product


def calculate_signatures(kmer_length,contigs):
    signatures = []
    kmers = possible_kmers(kmer_length)
    for c in contigs:
        signature = ss(kmer_length,c, kmers).kmer_frequencies
        signatures.append(signature)
    return signatures


def possible_kmers(k):
    kmer_dict = {}
    kmer_list = [''.join(x) for x in product('ATGC', repeat=k)]
    for i in range(len(kmer_list)):
        kmer_dict[kmer_list[i]] = i
    return kmer_dict
