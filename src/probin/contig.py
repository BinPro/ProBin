class Contig():
    def __init__(self,id=None,seq=None,kmer_len=4):
        self.id = id
        self.seq = seq.upper()
        self.kmer_len = kmer_len
        self.profile = self.generate_profile()
    def generate_profile(self):
        if len(self.seq) < self.kmer_len:
            raise Exception("Sequence too short")
        