import numpy as np

def read_mappings_to_log_coverage(rm, length,read_length):
    return np.log(0.1+rm.astype(float)*read_length/length)

def relative_frequency_to_log_coverage(rf,length,total_nu_bases):
    return np.log(0.1+rf*total_nu_bases/length)
