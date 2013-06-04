import os
from argparse import ArgumentParser
def main_parser():
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    subparsers = parser.add_subparsers()

    parser_bin = subparsers.add_parser('bin')
    parser_bin.add_argument('file', 
        help='specify input file on FASTA format')
    parser_bin.add_argument('-o', '--output', default=os.getcwd(),
        help='specify the output directory. The default is current directory')
    parser_bin.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser_bin.add_argument('-k', '--kmer', default=4, type=int,
        help='specify the length of kmer to use, default 4')
    parser_bin.add_argument('-mc', '--model_composition', default='multinomial', type=str, choices=['multinomial','dirichlet'],
        help='specify the composition model to use, default multinomial.')
    parser_bin.add_argument('--model_coverage', default='None', type=str, choices=['isotropic_gaussian','None'],
        help='specify the abundance model to use, default: None')
    parser_bin.add_argument('-cf','--coverage_file',
        help='specify input file containing coverage information')

    parser_bin.add_argument('--first_data',
                            help='specify the name of the first column containing sample-data in the coverage file')
    parser_bin.add_argument('--last_data',
                            help='specify the name of the last column containing sample-data in the coverage file')

    parser_bin.add_argument('--read_length', type=int,
                            help='Specify the length of the reads, to enable coverage calculations')

    parser_bin.add_argument('-a', '--algorithm', default='em', type=str, choices=['kmeans','em'],
        help='specify the clustering algorithm to use, default em.')
    parser_bin.add_argument('-c', '--cluster_count', default=10, type=int,
        help='specify the number of cluster to use')
    parser_bin.add_argument('-i', '--iterations', default=150, type=int,
        help='specify the maximum number of iterations to allow before halting clustering')
    parser_bin.add_argument('-e', '--epsilon', default=1E-3, type=float,
        help='specify the precision of the clustering as stop condition')
    parser_bin.add_argument('-r', '--repeat', default=16, type=int,
        help='specify the number of times to run clustering with different start conditions')

    parser_bin.set_defaults(script='probin')

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
