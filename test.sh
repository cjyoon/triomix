# This script is provided to test TrioMix on a 1000 genomes family. 
# simply runnin `sh test.sh` would execute download, indexing of the data, simulation of contaminations, and executing TrioMix
# requires samtools to be installed and executable in PATH. 

# download M008 family's WGS from 1000 genomes project ftp.
proband
curl ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989418/NA19662.final.cram  -o proband.cram
curl ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989418/NA19662.final.cram.crai  -o proband.cram.crai

# father
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239902/NA19661.final.cram -O father.cram
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239902/NA19661.final.cram.crai -O father.cram.crai

# mother
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989417/NA19660.final.cram -O mother.cram
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989417/NA19660.final.cram.crai -O mother.cram.crai

# sibling
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989425/NA19685.final.cram -O sibling.cram
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989425/NA19685.final.cram.crai -O sibling.cram.crai

# download the reference fasta file
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta



##################################
# sibling contamination in offspring simulation. Calculating the total read in the CRAM file can take a while.
python simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram \
	-r 0 0 0.75 0.25 -o sibling25_in_offspring 

# run TrioMix on sibling contaminated data
python triomix.py -f father.cram \
	-m mother.cram \
	-c sibling25_in_offspring/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 4 \
	-s common_snp/grch38_common_snp.bed.gz \
	-p sibling25_in_offspring -o results


##################################
# mother contamination by offspring simulation
python simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram \
	-r 0 0.25 0.75 0 -o mother25_in_offspring

# run TrioMix on mother contaminated data
python triomix.py -f father.cram \
	-m mother.cram \
	-c mother25_in_offspring/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 4 \
	-s common_snp/grch38_common_snp.bed.gz \
	-p mother25_in_offspring -o results



##################################
# father contamination in offspring simulation
python simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0.25 0 0.75 0 -o father25_in_offspring

# run TrioMix on mother contaminated data
python triomix.py -f father.cram \
	-m mother.cram \
	-c father25_in_offspring/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 4 \
	-s common_snp/grch38_common_snp.bed.gz \
	-p father25_in_offspring -o results


##################################
# mother contamination in father simulation
python simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0.75 0.25 0 0 -o mother25_in_father 

# run TrioMix on mother contaminated data
python triomix.py -f mother25_in_father/familymix.bam \
	-m mother.cram \
	-c proband.cram \
	-r Homo_sapiens_assembly38.fasta -t 4 \
	-s common_snp/grch38_common_snp.bed.gz \
	-p mother25_in_father -o results \
	--parent



##################################
# multiple contamination simulation, father=10%, mother=20%, offspring=40%, sibling=30%
python simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0.10 0.20 0.40 0.30 -o complexmix 

# run TrioMix on mother contaminated data
python triomix.py -f father.cram \
	-m mother.cram \
	-c complexmix/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 4 \
	-s common_snp/grch38_common_snp.bed.gz \
	-p complexmix -o results


##################################
# offspring contamination in father simulation
python simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0.75 0 0.25 0 -o offspring25_in_father

# run TrioMix on mother contaminated data
python triomix.py -f offspring25_in_father/familymix.bam \
	-m mother.cram \
	-c proband.cram \
	-r Homo_sapiens_assembly38.fasta -t 4 \
	-s common_snp/grch38_common_snp.bed.gz \
	-p offspring25_in_father -o results \
	--parent 

