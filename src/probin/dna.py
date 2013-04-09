from itertools import product
from collections import Counter, defaultdict
from random import randint
import sys
import os
import numpy as np

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

    def __init__(self,id,seq,phylo=None,calc_sign=False):
        if not self.kmer_len:
            raise Exception("Please run DNA.generate_kmer_hash(kmer_len) first.")
        self.id = id
        self.phylo = phylo
        self.seq = seq.upper().split("N")
        self.signature = None
        if calc_sign:
            self.calculate_signature()

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

    def split_seq_random(self,l,n):
        # left_overs contains tuples of left-over sequences 
        # of the genome and its starting position
        left_overs = [(self.full_seq,0)]
        parts = []
        i = 0
        left_over_lengths = [len(part[0]) for part in left_overs]
        max_length = left_over_lengths[0]
        while (i < n and max_length>l):
            # Calculate the lengths of all parts
            # Randomly choose which seq-part to sample
            index = np.argmax(left_over_lengths)
            
            if left_over_lengths[index]>l:
                part, new_left_overs = self._pick_part(
                    left_overs[index][0],left_overs[index][1],left_over_lengths[index],l)
                i+=1
                part_seq = DNA(id = self.id, seq=part[0])
                part_seq.start_position = part[1]
                parts.append(part_seq)
                left_overs += new_left_overs

                left_over_lengths = [len(part[0]) for part in left_overs]
                max_length = max(left_over_lengths)
        return parts

    def split_seq(self,l):
        # Parts contains tuples of sequences 
        # of the genome and its starting position
        seq = self.full_seq
        seq_l = len(seq)
        parts = []
        n = seq_l/l

        for i in xrange(n):
            pos = i*l
            part_seq = seq[pos:(pos+l)]
            part =  DNA(id = self.id, seq=part_seq)
            part.start_position = pos
            parts.append(part)
        part_seq = seq[n*l:]
        part = DNA(id =self.id, seq=part_seq)
        part.start_position = n*l
        parts.append(part)

        return parts

    def _pick_part(self, seq,seq_start,seq_l,l):
        start = randint(0, (seq_l-l))
        end = start+l
        part = (seq[start:end],seq_start + start)
        left = (seq[0:start],seq_start)
        right = (seq[end:],seq_start+end)
        return part, [left,right]
