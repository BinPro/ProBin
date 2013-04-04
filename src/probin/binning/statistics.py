# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:07:24 2013
Calculate recall and precision for clustering algorithms. Clustering 
@author: Brynjar Smári Bjarnason
"""
from pandas import DataFrame, Series
import sys

def recall(clustering, contigs):
    rows = []
    for contig in contigs:
        cluster = None
        for i,c in enumerate(clustering):
            if contig in c:
                if cluster != None:
                    raise Exception("Contig in multiple clusters, not good!")
                cluster = i
                break
        if cluster == None:
            raise Exception("Contig not in any cluster, not good!")
        serie = Series([contig.phylo_family,contig.phylo_genus,contig.phylo_species,cluster,len(contig.full_seq)],index=["family","genus","species","cluster","contig_size"],name=contig.id)
        rows.append(serie)
    df = DataFrame(rows)
    print>>sys.stderr,  df.head()
    seq_species = df.groupby([df.family,df.genus,df.species])
    recall = {}
    for name,group in seq_species:
        tot_nucleotides = group.sum()["contig_size"]
        recall[group.index[0]] = group.groupby(df.cluster).sum().max()["contig_size"] / tot_nucleotides
    return recall



def _print_confusion_matrix(clustering, matrix):
    pass