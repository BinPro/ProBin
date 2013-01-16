#!/usr/bin/env python
class Kmer_composition:
    def __init__(self,kmer_length, contig):
        self.kmer_length = kmer_length
        self.kmer_frequencies = self.init_freq(contig)
    def __len__(self):
        return self.kmer_length
    def frequencies(self):
        return self.kmer_frequencies
    def init_freq(self,contig):
        l = self.kmer_length
        sub = "A"*l
        freq = self.occurrences(contig.seq,sub)
        return [freq]
    def occurrences(self,string, sub):
        count = start = 0
        while True:
            start = string.find(sub, start) + 1
            if start > 0:
                count+=1
            else:
                return count
