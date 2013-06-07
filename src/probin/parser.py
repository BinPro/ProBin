import os
from argparse import ArgumentParser
def main_parser():
    #=============================
    #Default parameters for all bin subparsers
    #=============================   
    default_bin_parser = ArgumentParser(add_help=False)
    
    default_bin_parser.add_argument('-a', '--algorithm', default='em', type=str, choices=['kmeans','em'],
        help='specify the clustering algorithm to use, default em.')
    default_bin_parser.add_argument('-c', '--cluster_count', default=10, type=int,
        help='specify the number of cluster to use')
    default_bin_parser.add_argument('-r', '--runs', default=16, type=int,
        help='specify the number of times to run clustering with different start conditions')
    default_bin_parser.add_argument('-i', '--iterations', default=150, type=int,
        help='specify the maximum number of iterations in each run allowed before halting clustering')
    default_bin_parser.add_argument('-e', '--epsilon', default=1E-3, type=float,
        help='specify the log precision of the clustering as stop condition')
    default_bin_parser.add_argument('-o', '--output', default=os.getcwd(),
        help='specify the output directory. The default is current directory')
    default_bin_parser.add_argument('-v', '--verbose', action='store_true', default=False,
        help='information written to stderr during execution.')
    default_bin_parser.add_argument('-cent','--centroids', default=None,
        help='specify predefined centroids (NOT IMPLEMENTET YET)')
    default_bin_parser.add_argument('-s', '--serial', action='store_true',
        help='execute clustering serial, for debug')
    
    
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    subparsers = parser.add_subparsers(title="Select bin for clustering or preprocess for converting coverage to usable format")
    #=============================
    #Clustering parsers
    #=============================    
    parser_bin = subparsers.add_parser('bin')
    parser_bin.set_defaults(script='probin')
    #=============================
    #Select one type of model you will use
    #=============================    
    bin_subparsers = parser_bin.add_subparsers()
    #=============================
    #Composition parameters
    #=============================
    parser_composition = bin_subparsers.add_parser('composition',parents=[default_bin_parser])
    parser_composition.add_argument('file', 
        help='specify input file on FASTA format')
    parser_composition.add_argument('-m','--model', default="multinomial", type=str, choices=['multinomial','dirichlet'],
        help='specify the composition model to use, default multinomial.')
    parser_composition.add_argument('-k', '--kmer', default=4, type=int,
        help='specify the length of kmer to use, default 4')
    parser_composition.set_defaults(model_type='composition')

    #=============================
    #Coverage parameters
    #=============================
    parser_coverage = bin_subparsers.add_parser('coverage',parents=[default_bin_parser])
    parser_coverage.add_argument('file',
        help='specify input file containing coverage information')
    parser_coverage.add_argument('first_data',
        help='specify the name of the first column containing sample-data in the coverage file')
    parser_coverage.add_argument('last_data',
        help='specify the name of the last column containing sample-data in the coverage file')
    parser_coverage.add_argument('-m','--model', default="isotropic_gaussian", type=str, choices=['isotropic_gaussian'],
        help='specify the abundance model to use, default: isotropic_gaussian')

    parser_coverage.add_argument('--read_length', type=int, default=100,
        help='Specify the length of the reads, to enable coverage calculations, default 100')

    parser_coverage.set_defaults(model_type='coverage')


    #=============================
    #preprocessing steps
    #=============================
    parser_preprocess = subparsers.add_parser('preprocess')
    
    parser_preprocess.add_argument('files',nargs='*',
                                   help='speciy input files in masmvali output format')

    parser_preprocess.add_argument('-f','--format', 
                                   choices=['masmvali','bam'],
                                   default='bam',
                                   help='specify format of input files, default bam.')
    parser_preprocess.add_argument('-o', '--output', 
        help='specify the output file.  The default is stdout')
    parser_preprocess.add_argument('-c', '--contigs', 
        help='specify the file containing the contigs.')
    parser_preprocess.add_argument('-s', '--strain', action='store_true',
                                   help="""Use this option if the strain 
that the reads unambiguosly mapped to should be included in output""")
    
    parser_preprocess.set_defaults(script='preprocess')

    return parser
