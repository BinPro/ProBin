# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 19:08:24 2013

@author: binni
"""
import sys
import numpy as np
from multiprocessing import Pool


def cluster(cluster_func, contigs, p, K, epsilon, iterations, runs, verbose, serial, processors, **kwargs):
    (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)
    results=[]
    for k in K:
        params = [(cluster_func, contigs, p, k, epsilon, iterations, kwargs) for run in xrange(runs)]
        if not serial:
            try:
                pool = Pool(processes=processors)
                runs_results = pool.map(_clustering_wrapper, params)
            except Exception as e:
                print >> sys.stderr, "EM clustering failed. Error: {0}, message: {1}".format(e,e.message)
                sys.exit(-1)
            finally:
                pool.close()
        else:
            runs_results = [_clustering_wrapper(param) for param in params]

        results.append((k,max(runs_results,key=lambda x: x[1])))
    if kwargs["BIC"]:
        #Do bic calculations
        from probin.output import Output
        N,D = contigs.shape
        bics = [(k,-2*bic[2]+k*D*np.log(N)) for k,bic in results]
        Output.write_bic(bics)
    
    return results

def _clustering_wrapper(params):
    cluster_func = params[0]
    return cluster_func(*params[1:-1],**params[-1])