# 2019.07.30
# creating merged bam from two sibling bam files with specified subsample (mixing) ratios


import os

# bams = ['/home/users/cjyoon/Projects/rheum/bam/PDP1.bam', '/home/users/cjyoon/Projects/rheum/bam/PDP4.bam']
configfile: 'sample_config.yaml'
chimeric_ratio = [str(i) for i in range(0, 11)]
print(config['sample'])
rule all:
    input:
        expand('final/{sample}.{ratio}.bam', ratio=chimeric_ratio, sample=config['sample'])

rule subsample:
    input:
        sibling1 = lambda wildcards: config['sample'][wildcards.sample]['sibling1'],
        sibling2 = lambda wildcards: config['sample'][wildcards.sample]['sibling2'],
    output:
        sibling1_sub = 'subsampled/{sample}.sib1.{ratio}.bam', 
        sibling2_sub = 'subsampled/{sample}.sib2.{ratio}.bam', 
    params:
        ratio = '{ratio}', 
    threads: 4
    run:
        sibling_ratio = (10-int(params.ratio))/10
        proband_ratio = (int(params.ratio)/10)
        cmd = f'sambamba view -t {threads} -f bam -h -s {proband_ratio} {input.sibling1} > {output.sibling1_sub}; '
        print(cmd)
        os.system(cmd)
        cmd = f'sambamba view --subsampling-seed 1 -t {threads} -f bam -h -s {sibling_ratio} {input.sibling2} > {output.sibling2_sub}; '
        print(cmd)
        os.system(cmd)


rule merge:
    input:
        sibling1_sub = 'subsampled/{sample}.sib1.{ratio}.bam', 
        sibling2_sub = 'subsampled/{sample}.sib2.{ratio}.bam', 
    output:
        merged_bam = 'merged/{sample}.{ratio}.bam'
    log:
        'logs/{sample}.{ratio}.merge.log'
    threads: 4
    shell:
        '(sambamba merge -t {threads} {output.merged_bam} {input.sibling1_sub} {input.sibling2_sub}) &> {log}'

rule readgroup:
    input:
        merged_bam = 'merged/{sample}.{ratio}.bam'
    output:
        readgroup_bam = 'final/{sample}.{ratio}.bam'
    threads: 4
    params:
        ratio =   '{ratio}', 
        sample = '{sample}'
    shell:
        'java -jar  ~/tools/picard.jar AddOrReplaceReadGroups I={input.merged_bam} O={output.readgroup_bam} RGID={params.sample}{params.ratio} RGLB={params.sample}{params.ratio} RGPL=ILLUMINA RGPU={params.sample}{params.ratio} RGSM={params.sample}{params.ratio}; '
        'sambamba index {output.readgroup_bam}'
