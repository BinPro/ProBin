import os
from Bio import SeqIO

def main_preprocess(args):
    """ Extracts information from several files and outputs as a single table.

    The files submitted contains information about mapping of reads to contigs.
    The contigs are contained in the contigs file.

    """
    # create a dictionary that will contain all contig ids
    contig_ids = {}
    seqs = list(SeqIO.parse(args.contigs,"fasta"))
    for seq in seqs:
        pass
    args.output.write("This is the output " + os.linesep)
