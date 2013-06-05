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

Dependencies
-----------
Jinja2==2.6
Pygments==1.6
argparse==1.2.1
biopython==1.61
distribute==0.6.35
docutils==0.10
ipython==0.13.2
line-profiler==1.0b3
matplotlib==1.2.1
nose==1.3.0
numpy==1.7.1
openpyxl==1.6.2
pandas==0.11.0
python-dateutil==2.1
pytz==2013b
pyzmq==13.1.0
scipy==0.12.0
six==1.3.0
stevedore==0.8
tornado==3.0.1
virtualenv==1.9.1
virtualenv-clone==0.2.4
virtualenvwrapper==3.7
wsgiref==0.1.2



TODO
====
- [ ] Use single model attribute choice (merge mc and model-coverage)
- [ ] Have the input for multinomial as parameter rather than required (might only want coverage)
- [ ] Standardized format for the coverage input file and the first,last data parsing.
- [ ] Make the preprocessing generate the standardized coverage input file
