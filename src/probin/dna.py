from itertools import product
from collections import Counter, defaultdict
from random import randint
import sys
import os

# optimized sliding window function from
# http://stackoverflow.com/a/7636587
from itertools import tee, izip

def window(seq,n):
    els = tee(seq,n)
    for i,el in enumerate(els):
        for _ in xrange(i):
            next(el, None)
    return izip(*els)

class DNA(object):
    BASE_COMPLEMENT = {"A":"T","T":"A","G":"C","C":"G"}
    kmer_hash={}
    kmer_len = None

    def __init__(self,id,seq):
        if not self.kmer_len:
            raise Exception("Please run DNA.generate_kmer_hash(kmer_len) first.")
        self.id = id
        self.seq = seq.upper().split("N")
        self.signature = None

    def __repr__(self):
        return "{0}".format(self.id)        
    
    @classmethod
    def generate_kmer_hash(cls,kmer_len):
        if cls.kmer_hash:
            raise Exception("Already initialized, can't change during execution")
        cls.kmer_len = kmer_len
        counter = 0
        for kmer in product("ATGC",repeat=cls.kmer_len):
            kmer= ''.join(kmer)
            if kmer not in cls.kmer_hash:
                cls.kmer_hash[kmer] = counter
                rev_compl = ''.join([cls.BASE_COMPLEMENT[x] for x in reversed(kmer)])
                cls.kmer_hash[rev_compl] = counter
                counter += 1
        cls.kmer_hash_count = counter
    
    @property
    def full_seq(self):
        return "N".join(self.seq)
    
    def calculate_signature(self):
        not_in_hash = 0
        self.signature = Counter()
        for fragment in self.seq:
            if len(fragment) < self.kmer_len:
                continue
            (indexes,not_in_hash) = self._get_kmer_indexes(fragment)
            self.signature.update(indexes)
            if not_in_hash:
                sys.stderr.write("Sequence id: %s, skipped %i kmers that were not in dictionary%s" % (self.id,not_in_hash,os.linesep)) 
    

    def pseudo_count(self,index):
        return self.signature.get(index,0)+1

    @property
    def pseudo_counts(self):
        for i in xrange(self.kmer_hash_count):
            yield self.pseudo_count(i)

    def _get_kmer_indexes(self,seq):
        indexes = defaultdict(int)
        not_in_hash = 0

        # This call is the most time consuming 
        # Therefore it is wise to keep it as fast as possible
        kmers = Counter(window(seq,self.kmer_len))

        # Once the kmers are counted, it's safe to do the checking
        for kmer_tuple, count in kmers.iteritems():
            kmer = "".join(kmer_tuple)
            if kmer in self.kmer_hash:
                indexes[self.kmer_hash[kmer]] += count
            else:
                not_in_hash += 1
        return (indexes,not_in_hash)

    def split_seq_to_signatures(self,l,n):
        gen_l = len(self.full_seq)
        parts = []
        for i in xrange(n):
            start = randint(0,(gen_l-l))
            end = start+l
            part = DNA(id = "nonsense!", seq=self.full_seq[start:end])
            part.calculate_signature()
            parts.append(part)
        return [part.signature for part in parts]
