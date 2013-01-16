#!/usr/bin/env python
from ProBin.Helpers.sequence_signature import Sequence_signature as ss


def calculate_signatures(kmer_length,contigs):
    signatures = []
    for c in contigs:
        signature = ss(kmer_length,c).kmer_frequencies
        signatures.append(signature)
    return signatures
