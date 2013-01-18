#!/usr/bin/env bash
echo "Sampling contigs"
./contig_generation.py ~/Dropbox/Shared/Binni_Johannes/references.fa -v -n 100 --max_length 10000 --min_length 30 -o ../results/generated_contigs.fa
echo "Continue with probin"
../../../src/ProBin.py ../results/generated_contigs.fa -o ../results/log_prob_for_contigs.fa