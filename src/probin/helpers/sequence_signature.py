from  array import array

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
        
        for i in range(len(s) - (k - 1)):
            kmer = s[i:i+k]
            if "N" not in kmer:
                kmers[possible_kmers[str(kmer.upper())]] += 1

        return kmers
