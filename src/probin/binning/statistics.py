"""
Created on Thu Mar 28 10:07:24 2013
Calculate recall and precision for clustering algorithms. Clustering 
@author: Brynjar Smari Bjarnason
"""
from pandas import DataFrame, Series, pivot_table
import numpy as np
    
def confusion_matrix(contigs,clustering,levels=["family","genus","species"]):
    _get_phylo(contigs)
    df = _create_dataframe(contigs,clustering)
    cm = [pivot_table(df,rows=levels[:i+1],cols=["cluster"],aggfunc=np.sum) for i in xrange(len(levels))]

    return cm

def recall(contigs,clustering):
    confusion_matrixes = confusion_matrix(contigs,clustering)
    recall = [matrix.div(matrix.sum(axis=1),axis=0) for matrix in confusion_matrixes]
    return recall

def precision(contigs,clustering):
    confusion_matrixes = confusion_matrix(contigs,clustering)
    precision = [matrix.div(matrix.sum()) for matrix in confusion_matrixes]
    return precision
def _create_dataframe(contigs,clustering):
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
    return DataFrame(rows)

def _get_phylo(contigs):
    for contig in contigs:
        (contig.phylo_family, contig.phylo_genus, contig.phylo_species) = contig.phylo.split("|")