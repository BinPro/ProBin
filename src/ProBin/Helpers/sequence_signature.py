class Sequence_signature:
    def __init__(self,kmer_length, contig):
        self.kmer_length = kmer_length
        self.kmer_frequencies = self.init_freq(contig)
    def init_freq(self,contig):
        l = self.kmer_length
        sub = "A"*l
        freqs = self.kmer_composition(contig.seq,l)
        return freqs
    def possible_kmers(k):
        return [''.join(x) for x in product('ATGC', repeat=k)]
    def kmer_composition(s, k):
        kmers = {}
        
        for kmer in possible_kmers(k):
            kmers[kmer] = 0
            
        for i in range(len(s) - (k - 1)):
            kmer = s[i:i+k]
            kmers[kmer] += 1

        return kmers
