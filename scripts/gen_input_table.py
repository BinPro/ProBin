#!/usr/bin/env python
"""
@author: inodb
"""
import sys
import argparse

from Bio import SeqIO
from Bio.SeqUtils import GC


def get_gc_and_len_dict(fastafile):
    """Creates a dictionary with the fasta id as key and GC and length as keys
    for the inner dictionary."""
    out_dict = {}

    for rec in SeqIO.parse(fastafile, "fasta"):
        out_dict[rec.id] = {}
        out_dict[rec.id]["GC"] = GC(rec.seq)
        out_dict[rec.id]["length"] = len(rec.seq)

    return out_dict


def get_bedcov_dict(bedcoveragefile):
    """Uses the BEDTools genomeCoverageBed histogram output to determine mean
    coverage and percentage covered for each contig.
    
    Returns dict with fasta id as key and percentage covered and cov_mean as
    keys for the inner dictionary."""
    out_dict = {}

    for line in open(bedcoveragefile):
        cols = line.split()

        try:
            d = out_dict[cols[0]]
        except KeyError:
            d = {}
            out_dict[cols[0]] = d

        if int(cols[1]) == 0:
            d["percentage_covered"] = 100 - float(cols[4]) * 100.0
        else:
            d["cov_mean"] = d.get("cov_mean", 0) + int(cols[1]) * float(cols[4])

    return out_dict


def generate_input_table_for_probin(fastafile, bedcovfiles):
    """Writes the input table for Probin for stdout. See hackathon google
    docs."""
    fastad = get_gc_and_len_dict(fastafile)
    bedcovdicts = [ get_bedcov_dict(bcf) for bcf in bedcovfiles ]
    
    # Header
    sys.stdout.write("%s\t%s\t%s" % ( 'contig', 'length', 'GC' ))
    for bcf in bedcovfiles:
        sys.stdout.write("\tcov_mean_%s\tpercentage_covered_%s\n" % (bcf, bcf))

    # Content
    for acc in fastad:
        # fasta stats
        sys.stdout.write("%s\t%s\t%s"  %
            (
                acc,
                fastad[acc]['length'],
                fastad[acc]['GC']
            )
        )

        # bed coverage stats
        for bcd in bedcovdicts:
            try:
                # Print cov mean
                # If no 0 bases with 0 coverage, then all bases in the contig are covered
                sys.stdout.write("\t%f\t%f" % (bcd[acc]["cov_mean"], bcd[acc].get("percentage_covered", 100)))
            except KeyError:
                # No reads mapped to this contig
                sys.stdout.write("\t0\t0")
        sys.stdout.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafile", help="Contigs fasta file")
    parser.add_argument("bedcovfiles", nargs='+', help="genomeCoverageBed histogram output")
    args = parser.parse_args()
    generate_input_table_for_probin(args.fastafile, args.bedcovfiles)
