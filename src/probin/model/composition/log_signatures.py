import numpy as np
def signatures_to_log(contigs):
    pseudo_contigs = contigs +1
    return np.log(pseudo_contigs / pseudo_contigs.sum(axis=1,keepdims=True))

