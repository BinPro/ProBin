#!/usr/bin/env bash
# This experiment aim is to compare how different kmer-lengths affect 
# the genomic profile specificity for different taxonomic levels.
echo "Sampling contigs uniformly from reference genomes"
./contig_generation.py ~/Dropbox/Shared/Binni_Johannes/references.fa -v -n 100 --max_length 10000 --min_length 30 -o ../results/generated_contigs.fa
# Each contig should be scored against each genome. 
echo "Fitting the multinomial parameters to each genome"
../../../src/ProBin.py ../results/generated_contigs.fa -o ../results/log_prob_for_contigs.fa
echo "Scoring contigs against its proper genome"
echo "Scoring contigs against all other genomes"