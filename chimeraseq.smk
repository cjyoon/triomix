
import pysam
import cyvcf2
import re
import subprocess
import shlex
import sys

configfile: 'family_config.yaml'
configfile: srcdir('path_config.yaml')

SAMTOOLS = config['SAMTOOLS']
JAVA = config['JAVA']
VARSCAN = config['VARSCAN']
VT = config['VT']
GATK = config['GATK']
BGZIP = config['BGZIP']
TABIX = config['TABIX']
BCFTOOLS = config['BCFTOOLS']
GET_CANDIDATE_LOCI_COUNTS = srcdir('get_candidate_loci_summary.py')
CHIMERA_LIKELIHOOD = srcdir('chimera_likelihood.R')

def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile, reference_filename=reference_fasta)
    name = bam.header['RG'][0]['SM']
    return name



if config['assembly']=='GRCh38':
    chromosomes = ['chr'+str(i) for i in list(range(1, 23)) + ['X', 'Y']]
    REFERENCE = config['GRCh38_REFERENCE']
    DBSNP = config['GRCh38_DBSNP']

elif config['assembly']=='GRCh37':
    chromosomes = [str(i) for i in list(range(1, 23)) + ['X', 'Y']]
    REFERENCE = config['GRCh37_REFERENCE']
    DBSNP = config['GRCh37_DBSNP']
else:
    print(config['assembly'] + ' not supported, exiting...')
    sys.exit(1)


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name

def bgzip_tabix(vcf_file):
    '''compress and index vcf file'''
    indexed_vcf = vcf_file + '.gz'
    cmd = f'{BGZIP} -f {vcf_file}'
    compress = subprocess.Popen(shlex.split(cmd))
    compress.wait()
    cmd = f'{TABIX} -f -p vcf {indexed_vcf}'
    tabix = subprocess.Popen(shlex.split(cmd))
    tabix.wait()
    return indexed_vcf


rule all:
    input:
        expand('varscan/{family_id}.varscan.acgt.vtdcn.dbsnpa.vcf.gz', family_id = config['family']) + \
        expand('count_summary/{family_id}.counts', family_id = config['family']) + \
        expand('mle/{family_id}.counts.mle.pdf', family_id = config['family'])


rule write_sample_name:
    input: 
        father_bam = lambda wildcards: config['family'][wildcards.family_id]['father'],
        mother_bam = lambda wildcards: config['family'][wildcards.family_id]['mother'],
        child_bam = lambda wildcards: config['family'][wildcards.family_id]['child']
    output:
        sample_list = 'samplelist/{family_id}.samplelist',
    run:
        with open(output.sample_list, 'w') as f:
            for individual in input.father_bam +  input.mother_bam + input.child_bam:
                f.write(sampleNameBam(individual) + '\n')


rule mpileup2snv:
    input: 
        father_bam = lambda wildcards: config['family'][wildcards.family_id]['father'],
        mother_bam = lambda wildcards: config['family'][wildcards.family_id]['mother'],
        child_bam = lambda wildcards: config['family'][wildcards.family_id]['child'], 
        sample_list = 'samplelist/{family_id}.samplelist',

    params:
        family_id = '{family_id}',
        snvvcf = 'varscan/{family_id}_tmp/{family_id}.{chrom}.varscan.snp.vcf', 
        chrom = '{chrom}'

    output: 
        snvvcf = 'varscan/{family_id}_tmp/{family_id}.{chrom}.varscan.snp.vcf.gz',     
    threads: 1
    log:
        "logs/{family_id}.{chrom}.mpileup2snv.log"
    shell:
        "({SAMTOOLS} mpileup -B -Q 20 -q 20 -r {params.chrom} -f {REFERENCE} {input.father_bam} {input.mother_bam} {input.child_bam} | {JAVA} -jar {VARSCAN} mpileup2snp --min-coverage 10 --mean-reads2 2  --min-var-freq 0.01 --p-value 0.99 --output-vcf --strand-filter 1 --vcf-sample-list {input.sample_list} > {params.snvvcf}; {BGZIP} {params.snvvcf}; {TABIX} -p vcf {output.snvvcf} ) &> {log}"


rule mpileup2indel:
    input: 
        father_bam = lambda wildcards: config['family'][wildcards.family_id]['father'],
        mother_bam = lambda wildcards: config['family'][wildcards.family_id]['mother'],
        child_bam = lambda wildcards: config['family'][wildcards.family_id]['child'],
        sample_list = 'samplelist/{family_id}.samplelist',
    params:
        family_id = '{family_id}', 
        indelvcf = 'varscan/{family_id}_tmp/{family_id}.{chrom}.varscan.indel.vcf', 
        chrom = '{chrom}'
    output: 
        indelvcf = 'varscan/{family_id}_tmp/{family_id}.{chrom}.varscan.indel.vcf.gz',     
    threads: 1
    log:
        "logs/{family_id}.{chrom}.mpileup2indel.log"
    shell:
        "({SAMTOOLS} mpileup -B -Q 20 -q 20 -r {params.chrom} -f {REFERENCE} {input.father_bam} {input.mother_bam} {input.child_bam} | {JAVA} -jar {VARSCAN} mpileup2indel --min-coverage 10 --mean-reads2 2  --min-var-freq 0.01 --p-value 0.99 --output-vcf --strand-filter 1 --vcf-sample-list {input.sample_list} > {params.indelvcf}; {BGZIP} {params.indelvcf}; {TABIX} -p vcf {output.indelvcf} ) &> {log}"


rule merge_snp_indel:
    input:
        indelvcf = expand('varscan/{{family_id}}_tmp/{{family_id}}.{chrom}.varscan.indel.vcf.gz', chrom=chromosomes),
        snvvcf = expand('varscan/{{family_id}}_tmp/{{family_id}}.{chrom}.varscan.snp.vcf.gz', chrom=chromosomes)
    output:
        family_vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.vcf.gz'
    threads:1
    log:
        "logs/{family_id}.merge.log"
    shell:
        "({BCFTOOLS} concat {input.indelvcf} {input.snvvcf} -o {output.family_vcf} -O z -a ) &> {log} "

rule vt_decompoase_normalize:
    input:
        family_vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.acgt.vcf.gz'
    output:
        vt_vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.acgt.vtdcn.vcf.gz'
    params:
        vt_vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.acgt.vtdcn.vcf'
    threads: 1
    log:
        "logs/{family_id}.vt.log"
    shell:
        "({VT} decompose -s  {input.family_vcf} | {VT} normalize -n -r {REFERENCE} - >  {params.vt_vcf}; "
        "{BGZIP} {params.vt_vcf}; {TABIX} -p vcf {output.vt_vcf}) &> {log}"

rule remove_non_atgc:
    input:
        vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.vcf.gz'
    output:
        vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.acgt.vcf.gz'
    params:
        vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.acgt.vcf'
    threads: 1 
    run:
        vcf_handle = cyvcf2.VCF(input.vcf)
        output_handle = cyvcf2.Writer(params.vcf, vcf_handle)
        for variant in vcf_handle:
            if re.search(r'[ATGC]+', variant.REF) and re.search(r'[ATGC]+', variant.ALT[0]):
                output_handle.write_record(variant)

        vcf_handle.close()
        output_handle.close()

        bgzip_tabix(params.vcf)


rule annotate_dbsnp:
    input:
        vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.acgt.vtdcn.vcf.gz'
    output:
        vcf = 'varscan/{family_id}.varscan.acgt.vtdcn.dbsnpa.vcf.gz'
    params:
        vcf = 'varscan/{family_id}_tmp/{family_id}.varscan.acgt.vtdcn.dbsnpa.vcf'
    log:
        "logs/{family_id}.annotate_dbsnp.log"
    threads: 1
    shell:
        "({BCFTOOLS} annotate -a {DBSNP} -c ID {input.vcf} -O z -o {output.vcf}; tabix -p vcf {output.vcf}) &> {log} "

rule get_count_summary:
    input:
        vcf = 'varscan/{family_id}.varscan.acgt.vtdcn.dbsnpa.vcf.gz', 
        father_bam = lambda wildcards: config['family'][wildcards.family_id]['father'],
        mother_bam = lambda wildcards: config['family'][wildcards.family_id]['mother'],
        child_bam = lambda wildcards: config['family'][wildcards.family_id]['child'],
    output:
        counts = 'count_summary/{family_id}.counts'
    threads: 1
    params:
        family_id = '{family_id}'
    run:
        father_id = sampleNameBam(input.father_bam[0])
        mother_id = sampleNameBam(input.mother_bam[0])
        child_id = sampleNameBam(input.child_bam[0])
        shell("python {GET_CANDIDATE_LOCI_COUNTS} -i {input.vcf} -f {father_id} -m {mother_id} -c {child_id} -o count_summary -p {params.family_id}")

rule mle:
    input:
        counts = 'count_summary/{family_id}.counts'
    output:
        mle = 'mle/{family_id}.counts.mle.pdf'
    threads: 1
    shell:
        "Rscript {CHIMERA_LIKELIHOOD} -i {input.counts} -o mle"


