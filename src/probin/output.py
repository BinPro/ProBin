# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:10:45 2013

@author: binni
"""
import os
import sys
from datetime import datetime

class Output(object):
    file_name=None
    path=None
    start_time = None
    
    @classmethod
    def set_output_path(self,data_file,args):
        self.start_time = datetime.now()
        if not os.path.isdir(os.path.abspath(args.output)):
            self.path = os.getcwd()
        else:
            self.path = args.output
        self.file_name = "{0}_k{1}_c{2}_{3}".format(os.path.basename(data_file),args.kmer,args.cluster_count,args.algorithm)
        if os.path.isfile(os.path.join(self.path,self.file_name)):
            self.file_name = "{0}_no_override_{1}".format(self.file_name,datetime.now().strftime("%Y-%m-%d-%H.%M.%S"))
        print >> sys.stderr, "Result files created in {0} as {1}:".format(self.path,self.file_name)
    
    @classmethod
    def write_clustering_result(self,clusters, cluster_evaluation, centroids, arguments="", tmpfile=False,tmpfile_suffix=""):
        if not self.path or not self.file_name:
            print >> sys.stderr, "You need to call Output.set_output_path to initialize output"
            return None
        #CLUSTERING INFORMATION OUTPUT
        RESULT=["#Start time: {start_time}, now: {curr_time}, diff: {diff_time}",
                "#Clustering based on parameters: {args}",
                "#Result written to file: {directory}",
                "#Clustering evaluation: {clust_prob}",
                "#Cluster sizes",
                "{cluster_freq}",
                "#Output below is on the form",
                "#cluster_id,contig_n,contig_m,...",
                "#{divide}",
                "{clusters}"]
        curr_time = datetime.now()
        repr_centroids = ["#Centroid {0},{1}".format(i,",".join(map(str,centroid))) for i,centroid in enumerate(centroids)]
        cluster_sizes = [len(c) for c in clusters]
        tot_c = float(sum(cluster_sizes))
        cluster_freq = ["#Cluster {0}\t{1}\t{2}".format(i,c,c/tot_c) for i,c in enumerate(cluster_sizes)]
        cluster_contigs_id =  [ "{0},{1}".format(i,",".join([contig.id for contig in cluster]) )  for i,cluster in enumerate(clusters)]
        params =   {"args":arguments, "clust_prob":cluster_evaluation,
                    "centroids":os.linesep.join(repr_centroids),
                    "divide":"="*70,
                    "directory":self.file_name,
                    "cluster_freq":os.linesep.join(cluster_freq),
                    "clusters":os.linesep.join(cluster_contigs_id),
                    "start_time":self.start_time,
                    "curr_time":curr_time,
                    "diff_time":(curr_time-self.start_time)}
        if tmpfile:
            outfile = os.path.join(self.path,"_".join([self.file_name,str(tmpfile_suffix), curr_time.time().strftime("%H.%M.%S.%f")]))
        else:
            outfile = os.path.join(self.path,self.file_name)
        with open(outfile,"w") as clustinf:
            clustinf.write(os.linesep.join(RESULT).format(**params))
