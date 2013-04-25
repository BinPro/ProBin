from argparse import ArgumentParser
def main_parser():
    parser = ArgumentParser(description="Clustering of metagenomic contigs")
    subparsers = parser.add_subparsers()

    parser_bin = subparsers.add_parser('bin')
    parser_bin.add_argument('files', nargs='*', 
        help='specify input files on FASTA format, default is stdin')
    parser_bin.add_argument('-o', '--output', 
        help='specify the output file.  The default is stdout')
    parser_bin.add_argument('-v', '--verbose', action='store_true',
        help='information written to stderr during execution.')
    parser_bin.add_argument('-k', '--kmer', default=4, type=int,
        help='specify the length of kmer to use, default 4')
    parser_bin.add_argument('-mc', '--model_composition', default='multinomial', type=str, choices=['multinomial','dirichlet'],
        help='specify the composition model to use, default multinomial.')
    parser_bin.add_argument('-a', '--algorithm', default='em', type=str, choices=['kmeans','em'],
        help='specify the clustering algorithm to use, default em.')
    parser_bin.add_argument('-c', '--cluster_count', default=10, type=int,
        help='specify the number of cluster to use')

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
