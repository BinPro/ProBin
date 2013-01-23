from  array import array
from collections import defaultdict
from  itertools import product

class SequenceSignature:
    def __init__(self,kmer_length, contig, possible_kmers):
        self.kmer_length = kmer_length
        self.kmer_frequencies = self.init_freq(contig, possible_kmers)
    def init_freq(self,contig, possible_kmers):
        l = self.kmer_length
        freqs = self.kmer_composition(contig.seq,l, possible_kmers)
        return freqs
    def kmer_composition(self, s, k, possible_kmers):
        kmers = array('H', [0]*len(possible_kmers))
        s = s.upper()
        for i in range(len(s) - (k - 1)):
            kmer = s[i:i+k]
            kmers[possible_kmers[str(kmer)]] += 1
        return kmers

def possible_kmers(k):
    kmer_dict = defaultdict(int)
    kmer_list = [''.join(x) for x in product('ATGC', repeat=k)]
    kmer_list.insert(0,'REST')
    for i in range(len(kmer_list)):
        kmer_dict[kmer_list[i]] = i
    return kmer_dict


