# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 19:08:24 2013

@author: binni
"""
import sys
import numpy as np
from multiprocessing import Pool, cpu_count


def cluster(cluster_func, contigs, p, K, epsilon, iterations, runs, verbose, serial, **kwargs):
    (max_clusters, max_clustering_prob,max_centroids) = (None, -np.inf, None)
    params = [(cluster_func, contigs, p, K, epsilon, iterations, kwargs) for run in xrange(runs)]
    if not serial:
        try:
            pool = Pool(processes=cpu_count())
            results = pool.map(_clustering_wrapper, params)
        except Exception as e:
            print >> sys.stderr, "EM clustering failed. Error: {0}, message: {1}".format(e,e.message)
            sys.exit(-1)
        finally:
            pool.close()
    else:
        results = [_clustering_wrapper(param) for param in params]
        
    return max(results,key=lambda x: x[1])

def _clustering_wrapper(params):
    cluster_func = params[0]
    return cluster_func(*params[1:-1],**params[-1])