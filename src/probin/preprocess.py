import os
from Bio import SeqIO
import pandas as p

def main_preprocess(args):
    """ Extracts information from several files and outputs as a single table.

    The files submitted contains information about mapping of reads to contigs.
    The contigs are contained in the contigs file. Columns in resulting
    file will have names from the input files, make sure the names are unique.

    """
    # create a dictionary that will contain all contig ids
    contig_ids = []
    contig_lengths = []
    seqs = list(SeqIO.parse(args.contigs,"fasta"))
    for seq in seqs:
        contig_ids.append(seq.id)
        contig_lengths.append(len(seq.seq))
    series_dict = {}
    series_dict['contig_length'] = p.Series(contig_lengths,index=contig_ids)
    if args.format == 'masmvali':
        for input_file in args.files:
            file_base_name = ".".join(os.path.basename(input_file).split(".")[0:-1])
            input_df = p.io.parsers.read_table(input_file,
                                               sep='\t',
                                               index_col=0)
            series_dict[file_base_name] = input_df.unamb_tot_nr_reads
            if args.strain:
                series_dict[file_base_name+'_strain'] = input_df.unamb_dominant_strain
            if args.total_read:
                series_dict['total_read_count'] = input_df.tot_nr_reads
            if args.ratio:
                series_dict['contigs_genome_length_ratio'] = input_df.contigs_genome_length_ratio
        
        df = p.DataFrame(series_dict)
        df.to_csv(args.output, sep='\t')
