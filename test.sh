# This script is provided to test TrioMix on a 1000 genomes family. 
# simply runnin `sh test.sh` would execute download, indexing of the data, simulation of contaminations, and executing TrioMix
# requires SAMBAMBA to be installed and executable. 

# download M008 family's WGS from 1000 genomes project ftp.
# #proband
# curl ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989418/NA19662.final.cram  -o NA19662.final.cram
# #paternal 
# curl ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239902/NA19661.final.cram -o NA19661.final.cram
# #maternal
# curl ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989417/NA19660.final.cram -o NA19660.final.cram
# #sibling
# curl ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989425/NA19685.final.cram -o NA19685.final.cram

# #download the reference fasta file
# wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
# samtools faidx Homo_sapiens_assembly38.fasta

# #rename so that easier to work with
# mv NA19662.final.cram M008_proband.cram 
# mv NA19661.final.cram M008_father.cram
# mv NA19660.final.cram M008_mother.cram
# mv NA19685.final.cram M008_sibling.cram

# # convert cram to bam files


# samtools view -@ 4  -hbo M008_father.bam M008_father.cram
# samtools view -@ 4  -hbo M008_mother.bam M008_mother.cram
# samtools view -@ 4  -hbo M008_proband.bam M008_proband.cram
# samtools view -@ 4  -hbo M008_sibling.bam M008_sibling.cram

# # index the bam files
# for i in *.bam; do echo $i; samtools index -@ 4 $i; done

##################################
# sibling contamination simulation. Calculating the total read in the CRAM file can take a while.
# python simulate_familial_mixture.py \
# 	-f M008_father.bam \
# 	-m M008_mother.bam \
# 	-c M008_proband.bam \
# 	-s M008_sibling.bam \
# 	-r 0 0 0.75 0.25 -o M008_sibling25 

# run TrioMix on sibling contaminated data
# python triomix.py -f M008_father.bam \
# -m M008_mother.bam \
# -c M008_sibling25/familymix.bam \
# -r Homo_sapiens_assembly38.fasta -t 4 \
# -s common_snp/grch38_common_snp.bed.gz \
# -p sibling25 -o results


##################################
# mother contamination simulation
# python simulate_familial_mixture.py \
# -f M008_father.bam \
# -m M008_mother.bam \
# -c M008_proband.bam \
# -s M008_sibling.bam \
# -r 0 0.25 0.75 0 -o M008_mother25

# run TrioMix on mother contaminated data
# python triomix.py -f M008_father.bam \
# -m M008_mother.bam \
# -c M008_mother25/familymix.bam \
# -r Homo_sapiens_assembly38.fasta -t 4 \
# -s common_snp/grch38_common_snp.bed.gz \
# -p mother25 -o results



##################################
# father contamination simulation
# python simulate_familial_mixture.py \
# -f M008_father.bam \
# -m M008_mother.bam \
# -c M008_proband.bam \
# -s M008_sibling.bam  \
# -r 0.25 0 0.75 0 -o M008_father25 

# run TrioMix on mother contaminated data
# python triomix.py -f M008_father.bam \
# -m M008_mother.bam \
# -c M008_father25/familymix.bam \
# -r Homo_sapiens_assembly38.fasta -t 4 \
# -s common_snp/grch38_common_snp.bed.gz \
# -p father25 -o results


##################################
# multiple contamination simulation, father=10%, mother=20%, offspring=40%, sibling=30%
python simulate_familial_mixture.py \
-f M008_father.bam \
-m M008_mother.bam \
-c M008_proband.bam \
-s M008_sibling.bam  \
-r 0.10 0.20 0.40 0.30 -o M008_complexmix 

# run TrioMix on mother contaminated data
python triomix.py -f M008_father.bam \
-m M008_mother.bam \
-c M008_complexmix/familymix.bam \
-r Homo_sapiens_assembly38.fasta -t 4 \
-s common_snp/grch38_common_snp.bed.gz \
-p complexmix -o results



