# Triomix
Quantification of contamination or chimerism in whole-genome sequencing (WGS) data of parent-offspring trio. 


# Requirements
```bash
python v3.5 or later
- pysam

R v3.6.0 or later (including 4.0 or later)
- ggplot2
- tidyverse
- optparse
- bbmle
```
# Usage
`triomix.py` is the wrapper script that counts reference and variant read counts using `samtools mpileup` command, and then runs maximum-likelihood estimate to calculate the mixture. The tool can take a list of known common SNPs as BED file input, which can reduce the total run time. We have prepared a list of SNP sites for GRCh37 and GRCh38 reference sequence from `gnmoad v2` and `gnomad v3` respectively. These BED files are found in `common_snp` folder. Multithreading option is also supported with `-t` argument. 


```bash
# Whole-genome mode:
python triomix.py -f father.bam -m mother.bam -c child.bam -r reference.fasta -t 4

# Select snp mode:
python triomix.py -f father.bam -m mother.bam -c child.bam -r reference.fasta -t 4 -s common_snp/grch38_common_snp.bed.gz

```


```bash
$ python triomix.py -h
usage: triomix.py [-h] -f FATHER -m MOTHER -c CHILD -r REFERENCE [-s SNP] [-t THREAD] [-o OUTPUT_DIR] [-p PREFIX]
                  [--runmode {single,joint,all}] [-u {0,1}]

optional arguments:
  -h, --help            show this help message and exit
  -f FATHER, --father FATHER
                        Father's BAM or CRAM file
  -m MOTHER, --mother MOTHER
                        Mother's BAM or CRAM file
  -c CHILD, --child CHILD
                        Child's BAM or CRAM file
  -r REFERENCE, --reference REFERENCE
                        Reference FASTA file
  -s SNP, --snp SNP     Optional list of SNP sites as a BED (or BED.gz) file
  -t THREAD, --thread THREAD
                        Multithread to utilize. Default=1
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory. Default=current working directory
  -p PREFIX, --prefix PREFIX
                        prefix for the output file. If not specified, will use the SM tag from the child bam's header
  --runmode {single,joint,all}
                        Runmode for mle.R script. 'single' assumes only 1 contamination source within family. 'joint' calculates the fraction
                        of all family members jointly. 'all' runs both modes. Default=all
  -u {0,1}, --upd {0,1}
                        0: mle will filter out vaf=0 or 1 in sites where parental genotypes are homo-ref + homo-alt (GroupA SNPs) 1: mle will
                        identify UPDs which appears as contamination. Default=1

```
`triomix.py` internally calls `mle.R` and `plot_variants.R` to estimate DNA mixture and to plot the VAFs of variants.


# Output files
`*.counts` : raw table that contains the VAF of the parents and readcounts of the child at SNP positions. 

`*.counts.summary.tsv` : summary output of the `mle.R`. 

`*.counts.plot.pdf` : Depth and VAF plots of variants by `plot_variants.R`. 

`*.counts.*.homoalt.segments` : Segmentation of SNP VAFS in the child from GroupA SNPs. UPD will results in large segments with mean VAF value of 0 or 1. 

## SNP groups
Triomix classifies SNPs into three groups based on the parental genotypes. Each SNP types would have different patterns of VAFs in the offspring, which we used to infer the mixtures. 
```
GroupA: homo-ref + homo-alt (or vice versa) -> het (child)
GroupB: homo-ref + het (or vice versa) -> homo-ref or het (child)
GroupC: homo-ref + homo-ref -> homo-ref (child)
```
![SNP gruops](/images/snp_grouop.png "SNP groups")

## columns of *.counts.summary.tsv
```
input_file: Input *.counts file used for calculting DNA mixture
sibling_mix_j: Fraction of sibling's DNA mixture obtained from the 'joint' calculation mode. Calculated with GroupA and GroupB SNPs.
father_mix_j: Fraction of father's DNA mixture obtained from the 'joint' calculation mode. Calculated with GroupA and GroupB SNPs.
mother_mix_j: Fraction of mother's DNA mixture obtained from the 'joint' calculation mode. Calculated with GroupA and GroupB SNPs.
sibling_mix_s: Fraction of sibling's DNA mixture obtained from the 'single' calculation mode. Calculated with GroupB SNPs.
father_mix_s: Fraction of father's DNA mixture obtained from the 'single' calculation mode. Calculated with GroupA SNPs.
mother_mix_s: Fraction of mother's DNA mixture obtained from the 'single' calculation mode. Calculated with GroupA SNPs.
mendelian_error_rate: Fraction of alternative reads where both parents are homo-ref genotype. Calculated with GroupC SNPs.
```


# Docker
A Docker image is also available from Dockerhub `https://hub.docker.com/r/cjyoon/triomix/`. 
