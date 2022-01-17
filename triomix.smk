
import pysam
import cyvcf2
import re
import subprocess
import shlex
import sys

configfile: 'family_config.yaml'
configfile: srcdir('path_config.yaml')

TRIOMIX = srcdir('triomix.py')
TRIOMIX_PLOT = srcdir('plot_variant.R')

def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile, reference_filename=reference_fasta)
    name = bam.header['RG'][0]['SM']
    return name



if config['assembly']=='GRCh38':
    REFERENCE = config['GRCh38_REFERENCE']
    SNP_BED = '/home/users/cjyoon/scripts/chimeraseq/common_snp/grch38_common_snp.bed.gz'
elif config['assembly']=='GRCh37':
    REFERENCE = config['GRCh37_REFERENCE']
    SNP_BED = '/home/users/cjyoon/scripts/chimeraseq/common_snp/grch37_common_snp.bed.gz'
else:
    print(config['assembly'] + ' not supported, exiting...')
    sys.exit(1)


rule all:
    input:
        expand('triomix_snp/{family}.counts.summary.tsv', family = config['family']),
        expand('triomix_snp/{family}.counts.plot.pdf', family=config['family'])
        # expand('triomix_wgs/{family}.counts.summary.tsv', family = config['family'])  

rule triomix_snv:
    input: 
        father_bam = lambda wildcards: config['family'][wildcards.family]['father'],
        mother_bam = lambda wildcards: config['family'][wildcards.family]['mother'],
        child_bam = lambda wildcards: config['family'][wildcards.family]['child'], 

    params:
        family = '{family}',
        output_dir = 'triomix_snp', 
        snp_bed = SNP_BED

    output: 
        snvvcf = 'triomix_snp/{family}.counts.summary.tsv',  
        counts = 'triomix_snp/{family}.counts'   
    threads: 3
    log:
        "logs/{family}.triomix_snp.log"
    shell:
        "(python {TRIOMIX} -f {input.father_bam} -m {input.mother_bam} -c {input.child_bam} -o {params.output_dir} -s {params.snp_bed} -r {REFERENCE} -t {threads} -p {params.family} ) &> {log}"


rule triomix_vafplot:
    input: 
        counts = 'triomix_snp/{family}.counts',
        TRIOMIX_PLOT = TRIOMIX_PLOT
    params:
        family = '{family}',
        output_dir = 'triomix_snp', 
    output: 
        plot = 'triomix_snp/{family}.counts.plot.pdf',     
    threads: 1
    log:
        "logs/{family}.triomix.plot.log"
    shell:
        "(Rscript {input.TRIOMIX_PLOT} -i {input.counts} -o {params.output_dir} -r {REFERENCE} -s 0.1) &> {log}"




rule triomix_wgs:
    input: 
        father_bam = lambda wildcards: config['family'][wildcards.family]['father'],
        mother_bam = lambda wildcards: config['family'][wildcards.family]['mother'],
        child_bam = lambda wildcards: config['family'][wildcards.family]['child'], 

    params:
        family = '{family}',
        output_dir = 'triomix_wgs', 

    output: 
        snvvcf = 'triomix_wgs/{family}.counts.summary.tsv',     
    threads: 1
    log:
        "logs/{family}.triomix_wgs.log"
    shell:
        "(python {TRIOMIX} -f {input.father_bam} -m {input.mother_bam} -c {input.child_bam} -o {params.output_dir} -r {REFERENCE} -t {threads} -p {params.family} ) &> {log}"

