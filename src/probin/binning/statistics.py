# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:07:24 2013
Calculate recall and precision for clustering algorithms. Clustering 
@author: Brynjar Smári Bjarnason
"""
from pandas import DataFrame, Series
import sys

def recall(contigs,clustering):
    _get_phylo(contigs)
    df = _create_dataframe(contigs,clustering)

    phylo_calc = [df.groupby([df.family]),df.groupby([df.family,df.genus]),df.groupby([df.family,df.genus,df.species])]
    recall = {k:v for phylo in phylo_calc for k,v in _calc_recall(phylo).iteritems()}
    return recall

def precision(contigs,clustering):
    _get_phylo(contigs)
    df = _create_dataframe(contigs,clustering)
    
    cluster_groups = df.groupby(df.cluster)

    precision = {k:v for cluster_group in cluster_groups for k,v in _calc_precision(cluster_group).iteritems() }
    return precision
    
def confusion_matrix(contigs,clustering):
    _get_phylo(contigs)
    df = _create_dataframe(contigs,clustering)
    
    keys_family = df.groupby(df.family).groups.keys()
    data = [Series([0]*len(df.cluster.unique()),index=sorted(df.cluster.unique()),name=k) for k in keys_family]
    ndf = DataFrame(data)
    
    return None
    

def _calc_precision(cluster_group):
    name,cluster  = cluster_group
    cluster_size = cluster.contig_size.sum()
    precision = { \
                (name,cluster.family.name):\
                    cluster.groupby([cluster.family]).contig_size.sum().max() / cluster_size, \
                (name,cluster.family.name,cluster.genus.name):\
                    cluster.groupby([cluster.family,cluster.genus]).contig_size.sum().max() / cluster_size, \
                (name,cluster.family.name,cluster.genus.name,cluster.species.name):\
                    cluster.groupby([cluster.family,cluster.genus,cluster.species]).contig_size.sum().max() / cluster_size \
                }
    return precision

def _calc_recall(phylo_group_level):
    recall = {}
    for name,phylo_group in phylo_group_level:
        if type(name) is not tuple:
            name = tuple([name])
        tot_nucleotides = phylo_group.contig_size.sum()
        recall[name] = phylo_group.groupby(phylo_group.cluster).contig_size.sum().max() / tot_nucleotides
    return recall

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