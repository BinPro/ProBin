from  itertools import product

class SequenceSignature:
    def __init__(self,kmer_length, contig):
        self.kmer_length = kmer_length
        self.kmer_frequencies = self.init_freq(contig)
    def init_freq(self,contig):
        l = self.kmer_length
        sub = "A"*l
        freqs = self.kmer_composition(contig.seq,l)
        return freqs
    def possible_kmers(self, k):
        return [''.join(x) for x in product('ATGC', repeat=k)]
    def kmer_composition(self, s, k):
        kmers = {}
        
        for kmer in self.possible_kmers(k):
            kmers[kmer] = 0

        for i in range(len(s) - (k - 1)):
            kmer = s[i:i+k]
            if "N" not in kmer:
                kmers[str(kmer.upper())] += 1

        return kmers
