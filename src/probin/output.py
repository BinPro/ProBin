# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:10:45 2013

@author: binni
"""
import os
import sys
from datetime import datetime

class Output(object):
    path=None
    full_file_name=None
    start_time = None
    
    @classmethod
    def set_output_path(self,data_file,args):
        self.start_time = datetime.now()
        if not os.path.isdir(os.path.abspath(args.output)):
            self.path = os.getcwd()
        else:
            self.path = args.output
        self.full_file_name = "{0}_{1}".format(os.path.join(self.path,os.path.basename(data_file)),args.algorithm)
        if os.path.isfile(self.full_file_name):
            self.full_file_name = "{0}_no_override_{1}".format(self.full_file_name,datetime.now().strftime("%Y-%m-%d-%H.%M.%S"))
        print >> sys.stderr, "Result files created based on {0}:".format(self.full_file_name)
        
    @classmethod
    def write_clustering_result(self,clusters, cluster_evaluation, centroids, idx=None, arguments="", tmpfile=False,tmpfile_suffix="",cluster_number=None):
        if type(centroids) is tuple:
            centroids , sigma = centroids
        if not self.path or not self.full_file_name:
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
        cluster_sizes = [len(idx[clusters==i]) for i,c in enumerate(centroids)]
        tot_c = float(sum(cluster_sizes))
        cluster_freq = ["#Cluster {0}\t{1}\t{2}".format(i,c,c/tot_c) for i,c in enumerate(cluster_sizes)]
        cluster_contigs_id =  [ "{0},{1}".format(i,",".join(idx[clusters==i]) )  for i,c in enumerate(centroids)]
        params =   {"args":arguments, "clust_prob":cluster_evaluation,
                    "centroids":os.linesep.join(repr_centroids),
                    "divide":"="*70,
                    "directory":self.full_file_name,
                    "cluster_freq":os.linesep.join(cluster_freq),
                    "clusters":os.linesep.join(cluster_contigs_id),
                    "start_time":self.start_time,
                    "curr_time":curr_time,
                    "diff_time":(curr_time-self.start_time)}
        if tmpfile:
            outfile = "_".join([self.full_file_name,"tmp",str(tmpfile_suffix), curr_time.time().strftime("%H.%M.%S.%f")])
        else:
            outfile = "{0}_{1}".format(self.full_file_name,cluster_number)
        with open(outfile,"w") as clustinf:
            clustinf.write(os.linesep.join(RESULT).format(**params))
    @classmethod
    def write_bic(self,bics):
        outfile = "{0}_bic".format(self.full_file_name)
        with open(outfile,"w") as fh:
            fh.write("k,bic,-2*sum(bic),k*D*log(N),reg_prob\n")
            fh.writelines(["{0},{1},{2},{3},{4}{5}".format(k,bic,neg_bic,penalize,reg_prob,os.linesep) for k,bic,neg_bic,penalize,reg_prob in bics])
                
            