from itertools import product
from collections import Counter
class DNA(object):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash={}
    kmer_len = None
    @classmethod
    def generate_kmer_hash(cls,kmer_len):
        if cls.kmer_hash:
            raise Exception("Already initialized, can't change during execution")
        cls.kmer_len = kmer_len
        counter = 0
        for kmer in product("ATGC",repeat=cls.kmer_len):
            if kmer not in cls.kmer_hash:
                kmer = ''.join(kmer)
                cls.kmer_hash[kmer] = counter
                rev_compl = ''.join([cls.BASE_COMPLEMENT[x] for x in reversed(kmer)])
                cls.kmer_hash[rev_compl] = counter
                counter += 1
    
    def __init__(self,id,seq):
        if not self.kmer_len:
            raise Exception("Please run DNA.generate_kmer_hash(kmer_len) first.")
        self.id = id
        self.seq = seq.upper().split("N")
        self.signature = self.calculate_signature()
        
    def calculate_signature(self):
        signature = Counter()
        for fragment in self.seq:
            if len(fragment) < self.kmer_len:
                continue
            indexes = [self.kmer_hash[fragment[i:i+self.kmer_len]] for i in xrange(len(fragment) - (self.kmer_len-1))]
            signature.update(indexes)
        return signature

if __name__=="__main__":
    DNA.generate_kmer_hash(4)
    a = DNA(id="ADF",seq="ACTTNACTT")
    print a.signature
    b = DNA(id="ADFA",seq="ACTTTAAACCCACACACAACATTTGGAAAGGAGAGAGCCATTA")
    print b.signature
    c = DNA(id="ADADAD",seq='AAAATTTTACGTAGAGCCATTGAGACCTT')
    print c.signature
