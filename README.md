# ChimeraSeq
Identification of chimerism due to vanishing twin syndrome in whole-genome sequencing data of parent-offspring trio. 


# Requirements
```bash
Varscan v2.3.9 (https://github.com/dkoboldt/varscan)

vt (https://genome.sph.umich.edu/wiki/Vt)
python 3.5 or later
python packages
-cyvcf2
-pysam

R
R packages
-ggplot2
-tidyverse
-optparse
```
# One line command that execultes all using `Snakemake`
Use the supplied Snakefile `varscan_parallel.smk` to create one VCF per family. 
```bash
# If using GRCh38
snakemake -p -s chimeraseq.smk --config assembly=GRCh38
# If using GRCh37
snakemake -p -s chimeraseq.smk --config assembly=GRCh37
```
# Step by step guide
## Step1: Preparing input VCF file using Varscan2

If you prefer to use bash script to generate the VCF file, use the following commands. `$SAMPLE_LIST` file can be created by identifying `@RG SM:` tag for each bam file, and writing one line for each bam in the same order as used for `mpileup` input

```bash
# Prepare $SAMPLE_LIST file for Varscan
for bam in $FATHER_BAM $MOTHER_BAM $CHILD_BAM; do samtools view -H $bam  | grep '^@RG' | sed 's/.*SM://' | awk -F '\t' '{print $1}' >> $SAMPLE_LIST; done

# Varscan2 for SNVs
samtools mpileup -B -Q 20 -q 20 -f $REFERENCE $FATHER_BAM $MOTHER_BAM $CHILD_BAM | java -jar VarScan.v2.3.9.jar mpileup2snp --min-coverage 10 --mean-reads2 2  --min-var-freq 0.01 --p-value 0.99 --output-vcf --strand-filter 1 --vcf-sample-list $SAMPLE_LIST > family.snv.vcf; bgzip family.snv.vcf; tabix -p vcf family.snv.vcf.gz

# Varscan2 for indelx
samtools mpileup -B -Q 20 -q 20 -f $REFERENCE $FATHER_BAM $MOTHER_BAM $CHILD_BAM | java -jar VarScan.v2.3.9.jar mpileup2indel --min-coverage 10 --mean-reads2 2  --min-var-freq 0.01 --p-value 0.99 --output-vcf --strand-filter 1 --vcf-sample-list $SAMPLE_LIST > family.indel.vcf; bgzip family.indel.vcf; tabix -p vcf family.indel.vcf.gz

# combine indel and SNV vcf into a single vcf file
bcftools concat family.snv.vcf.gz family.indel.vcf.gz -o family.vcf.gz -O z -a 

# remove non ATGC bases from VCF file
python remove_nonacgt.py family.vcf.gz 

# Decompose and Normalize variants
vt decompose -s family.acgt.vcf.gz | vt normalize -r $REFERENCE - > family.acgt.vtdcn.vcf; bgzip family.acgt.vtdcn.vcf; tabix -p vcf family.acgt.vtdcn.vcf.gz

# annotate variants with dbSNP IDs
bcftools annotate -a $DBSNP -c ID family.acgt.vtdcn.vcf.gz -O z -o family.acgt.vtdcn.dbsnpa.vcf.gz

```

# Step2: Identify informative loci from Varscan VCF, and prepare a table for MLE estimates
From the VCF file, we can identify loci that is heterozygous in one of the parent and homozygous reference in the other parent. For each of these high-confidence loci, raw read depth and alternative read counts will be generated for downstream maximum likelihood estimate calculation. 

```bash
$ OUTPUT_DIR=result
$ PREFIX=family
$ python get_candidate_loci_summary.py -f father_rg -m mother_rg -c child_rg -o $OUTPUT_DIR -p $PREFIX
```
results in `{OUTPUT_DIR}/{PREFIX}.counts`. This is file is the sole input for the last step. 


# Step3: MLE estimate
```bash
$ Rscript chimera_likelihood.R -i $OUTPUT_DIR/$PREFIX.counts -o mle
```
Results in VAF histograms and MLE estimated results. 



