
import os
import re
SAMPLES = [re.sub(r'.vafs', '', i) for i in os.listdir() if i.endswith('vafs')]

rule all:
	input:
		expand('{sample}.vafs.hist.pdf', sample=SAMPLES) + 
		expand('{sample}.vafs.mle.pdf', sample=SAMPLES)


rule chimera_likelihood:
	input:
		inputfile = '{sample}.vafs'
	output:
		hist = '{sample}.vafs.hist.pdf', 
		mle = '{sample}.vafs.mle.pdf'
	shell:
		'Rscript /home/users/cjyoon/Projects/chimera/analysis/autism_trio/chimera_likelihood_20190807.R -i {input.inputfile}'