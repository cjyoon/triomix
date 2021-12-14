# Triomix
Quantification of contamination or chimerism in whole-genome sequencing (WGS) data of parent-offspring trio. 


# Requirements
```bash
python 3.5 or later
python packages
-pysam

R
R packages
-ggplot2
-tidyverse
-optparse
```
# Usage
`triomix.py` is the wrapper script that counts reference and variant read counts using `Varscan2`, normalizes the variant calls using `vt` and then runs maximum-likelihood estimate to calculate the mixture. The tool can also take list of known common SNPs as BED file input, which can reduce the total run time. We have prepared a list of SNP sites for GRCh37 and GRCh38 reference sequence from `gnmoad v2` and `gnomad v3` respectively. These BED files are found in `common_snp` folder. 


```bash
# Whole-genome mode:
python triomix.py -f father.bam -m mother.bam -c child.bam -r reference.fasta -t 4

# Select snp mode:
python triomix.py -f father.bam -m mother.bam -c child.bam -r reference.fasta -t 4 -s common_snp/grch38_common_snp.bed.gz

```


```bash
$ python triomix.py -h
usage: triomix.py [-h] -f FATHER -m MOTHER -c CHILD -r REFERENCE [-s SNP]
                  [-t THREAD] [-o OUTPUT_DIR] [-p PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -f FATHER, --father FATHER
                        Father's bam file
  -m MOTHER, --mother MOTHER
                        Mother's bam file
  -c CHILD, --child CHILD
                        Child's bam file
  -r REFERENCE, --reference REFERENCE
                        Reference fasta file
  -s SNP, --snp SNP     Optional list of SNP sites as a BED file
  -t THREAD, --thread THREAD
                        Multithread to utilize
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory
  -p PREFIX, --prefix PREFIX
                        prefix for the output file. If not specified, will use
                        the SM tag from the child's bam

```

# Docker
A Docker image is also available from Dockerhub `https://hub.docker.com/r/cjyoon/triomix/`. 
