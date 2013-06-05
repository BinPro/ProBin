ProBin
======

A program for binning metagenomic contigs to taxonomic rank by using nucleotide composition and correlation between samples. 


Execute ProBin
-------
```
- contigs.fna     contains the fasta formatted contigs
- k               the kmer size
- c               number of clusters
- r               the number of runs to execute clustering
- i               the number of maximum iterations per run of clustering
- e               the stop condition for the log differences between iterations in clustering
- a               algorithm to use (em,kmeans)
- mc              what model to use for composition
- cf              file with the coverage data
- model_coverage  what model to use for coverage
- first_data      the header of the first column in the coverage file
- last_data       the header of the last column in the coverage file
- read_length     the length of the reads for coverage calculations
- o               The directory where result files should be stored, otherwise current dir used
```

Example of executing ProBin
```
ProBin.py bin contigs.fna -mc multinomial -k 4 -c 10 -a em -r 10 -i 100 -e 0.001 \
              -cf coverage.tsv --model_coverage isotropic_gaussian --first_data 2012-03-25 \
              --last-data 2013-01-18 --read_length 100 -o /tmp/results

```




TODO
====
- [ ] Use single model attribute choice (merge mc and model-coverage)
- [ ] Have the input for multinomial as parameter rather than required (might only want coverage)
- [ ] Standardized format for the coverage input file and the first,last data parsing.
- [ ] Make the preprocessing generate the standardized coverage input file
