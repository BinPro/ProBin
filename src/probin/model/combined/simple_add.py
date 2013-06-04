import sys
from probin.model.composition import multinomial
from probin.model.coverage import isotropic_gaussian

def log_probability(seq,cov_matrix,prob_vector,mu,sigma):
    ig_p = isotropic_gaussian.log_pdf(cov_matrix,mu,sigma)
    mu_p = multinomial.log_probability(seq,prob_vector)
    return ig_p+mu_p


def fit_nonzero_parameters(dna_l,cov_matrix=None,expected_clustering=None,**kwargs):
    if not cov_matrix is None:
        if len(dna_l) != cov_matrix.shape[0]:
            sys.stderr.write("ERROR! Different numbers of contigs in fit nonzero parameters in simple add model!\n")
            sys.exit(-1)
        par_ig =  isotropic_gaussian.fit_nonzero_parameters(
            cov_matrix,
            expected_clustering=expected_clustering)
        par_mul = multinomial.fit_nonzero_parameters(
            dna_l,
            expected_clustering=expected_clustering)
    else:
        par_ig = (None,None)
        par_mul = multinomial.fit_nonzero_parameters(
            dna_l,
            expected_clustering=expected_clustering)

    par = (par_mul,par_ig[0],par_ig[1])

    return par

