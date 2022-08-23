# This script is provided to test TrioMix on a 1000 genomes family. 
# simply runnin `sh test.sh` would execute download, indexing of the data, simulation of contaminations, and executing TrioMix
# requires samtools to be installed and executable in PATH. 

# download M008 family's WGS from 1000 genomes project ftp.
# proband
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989418/NA19662.final.cram  -O proband.cram
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989418/NA19662.final.cram.crai  -O proband.cram.crai

# father
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239902/NA19661.final.cram -O father.cram
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239902/NA19661.final.cram.crai -O father.cram.crai

# mother
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989417/NA19660.final.cram -O mother.cram
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989417/NA19660.final.cram.crai -O mother.cram.crai

# sibling
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989425/NA19685.final.cram -O sibling.cram
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3989425/NA19685.final.cram.crai -O sibling.cram.crai

# download the reference fasta file
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
samtools faidx Homo_sapiens_assembly38.fasta


#################################
# get the location of github folder
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"


##################################
# offspring contaminated by a sibling simulation
python $SCRIPTPATH/simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram \
	-r 0 0 0.75 0.25 -o offspring75_sibling25

# run TrioMix on offspring contaminated by a sibling
python $SCRIPTPATH/triomix.py \
	-f father.cram \
	-m mother.cram \
	-c offspring75_sibling25/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 8 \
	-s $SCRIPTPATH/common_snp/grch38_common_snp.bed.gz \
	-p offspring75_sibling25 -o results


##################################
# offspring contaminated by mother simulation
python $SCRIPTPATH/simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram \
	-r 0 0.25 0.75 0 -o offspring75_mother25

# run TrioMix on offspring contaminated by mother
python $SCRIPTPATH/triomix.py \
	-f father.cram \
	-m mother.cram \
	-c offspring75_mother25/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 8 \
	-s $SCRIPTPATH/common_snp/grch38_common_snp.bed.gz \
	-p offspring75_mother25 -o results



##################################
# offspring contaminated by father simulation
python $SCRIPTPATH/simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0.25 0 0.75 0 -o offspring75_father25

# run TrioMix on offspring contaminated by father
python $SCRIPTPATH/triomix.py \
	-f father.cram \
	-m mother.cram \
	-c offspring75_father25/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 8 \
	-s $SCRIPTPATH/common_snp/grch38_common_snp.bed.gz \
	-p offspring75_father25 -o results



##################################
# multiple contamination simulation, father=10%, mother=20%, offspring=40%, sibling=30%
python $SCRIPTPATH/simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0.10 0.20 0.40 0.30 -o complexmix 

# run TrioMix on the complex contaminated case
python $SCRIPTPATH/triomix.py \
	-f father.cram \
	-m mother.cram \
	-c complexmix/familymix.bam \
	-r Homo_sapiens_assembly38.fasta -t 8 \
	-s $SCRIPTPATH/common_snp/grch38_common_snp.bed.gz \
	-p complexmix -o results


##################################
# mother contaminated by offspring simulation
python $SCRIPTPATH/simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0 0.75 0.25 0 -o mother75_offspring25

# run TrioMix on mother contaminated by offspring, parent mode
python $SCRIPTPATH/triomix.py \
	-f father.cram \
	-m mother75_offspring25/familymix.bam \
	-c proband.cram \
	-r Homo_sapiens_assembly38.fasta -t 8 \
	-s $SCRIPTPATH/common_snp/grch38_common_snp.bed.gz \
	-p mother75_offspring25 -o results \
	--parent 


##################################
# mother contaminated by father simulation
python $SCRIPTPATH/simulate_familial_mixture.py \
	-f father.cram \
	-m mother.cram \
	-c proband.cram \
	-s sibling.cram  \
	-r 0.25 0.75 0 0 -o mother75_father25 

# run TrioMix on mother contaminated by father, parent mode
python $SCRIPTPATH/triomix.py \
	-f father.cram \
	-m mother75_father25/familymix.bam  \
	-c proband.cram \
	-r Homo_sapiens_assembly38.fasta -t 8 \
	-s $SCRIPTPATH/common_snp/grch38_common_snp.bed.gz \
	-p mother75_father25 -o results \
	--parent

