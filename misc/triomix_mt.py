import sys
import os
import subprocess
import shlex
import re
import argparse
import json
import multiprocessing as mp
import pysam
import gzip

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--father', required=True, help="Father's bam file")
    parser.add_argument('-m', '--mother', required=True, help="Mother's bam file")
    parser.add_argument('-c', '--child', required=True, help="Child's bam file")
    parser.add_argument('-r', '--reference', required=True, help="Reference fasta file")
    parser.add_argument('--mtdna', required=True, help="contig  name of mitochondira in the reference genome")
    parser.add_argument('-t', '--thread', required=False, default=1, type=int, help="Multithread to utilize")
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-p', '--prefix', required=False, default=None, help="prefix for the output file. If not specified, will use the SM tag from the child's bam")
    args = vars(parser.parse_args())

    if args['prefix'] == None:
        args['prefix'] = sampleNameBam(args['child'])
    return args['father'], args['mother'], args['child'], args['reference'], args['mtdna'], args['thread'], args['output_dir'], args['prefix']


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name


def mpileup(father_bam, mother_bam, child_bam, region, output_dir):
    """Run mpileup in a defined region of interest"""
    
    child_id = sampleNameBam(child_bam)
    region_string = re.sub(r':|-', '_', region)
    # print(region_string)

    # # if SNP-BED file is provided 
    # if snp_bed != None:
    #     snp_bed_string = f' -l {snp_bed} '
    # else:
    #     snp_bed_string = ''

    # output_file = os.path.join(output_dir, f'{child_id}_{region_string}.mpileup')
    output_file_compressed = os.path.join(output_dir, f'{child_id}_{region_string}.mpileup.gz')
    print(output_file_compressed)

    mtdna_option = f'--max-depth 99999999' #added for mtdna
    if not os.path.isfile(output_file_compressed):
        cmd = f'{SAMTOOLS} mpileup {mtdna_option} -B -Q 20 -q 20 -r {region} -f {REFERENCE} {father_bam} {mother_bam} {child_bam} | gzip -f > {output_file_compressed}'
        os.system(cmd)
    # gzip_cmd = f'{GZIP} -f ' 
    # gzip_execute = subprocess.Popen(shlex.split(gzip_cmd), stdin=subprocess.PIPE, shell=True, stderr=subprocess.DEVNULL)
    
    
    # mpileup_cmd = f'{SAMTOOLS} mpileup -B -Q 20 -q 20 {snp_bed_string}-r {region} -f {REFERENCE} {father_bam} {mother_bam} {child_bam} '
    # mpileup_execute = subprocess.Popen(shlex.split(mpileup_cmd), stdout=gzip_execute.stdin, stderr=subprocess.DEVNULL, shell=True)
    
    # mpileup_execute.wait()
    # gzip_execute.stdin.close()
    # gzip_execute.wait()

    return  output_file_compressed



def vaf(alt_count, depth):
    if depth > 0:
        return float(alt_count/depth)
    else:
        return 0


def count_int(count):
    """handle '.' in alt/ref counts"""
    if count == None or count < 0:
        return 0
    else:
        return count

def parse_mpileup_mtdna(mpileup_line):
    """Parse mpileup result into counts string, modified for mtdna"""
    split_mpileup = mpileup_line.split('\t'); #print(split_mpileup)
    chrom = split_mpileup[0]
    pos = float(split_mpileup[1])
    ref = split_mpileup[2]

    totalDepth = 0 
    mismatches = 0
    totalCharacters = 0
    trio_alt_counts = dict({'father': None, 'mother': None, 'child': None})
    trio_depth_counts = dict({'father': None, 'mother': None, 'child': None})
    for individual, i in zip(['father', 'mother', 'child'], range(1, 4)):
        # initialize mismatch dict for each bam
        mismatch_dict = dict({'A': 0, 'C': 0, 'G': 0, 'T': 0, 'ins': 0, 'del': 0, 'depth': 0})

        bases_index = 3*i + 1
        depths_index = 3*i
        depths = int(split_mpileup[depths_index])
        totalDepth += depths
        mpiledup = split_mpileup[bases_index].upper()
        insertions = re.findall(r'\+[0-9]+[ACGTNacgtn]+', mpiledup)
        deletions = re.findall(r'-[0-9]+[ACGTNacgtn]+', mpiledup)

        mismatch_dict['ins'] = len(insertions)
        mismatch_dict['del'] = len(deletions)

        mpileupsnv = re.sub(r'\+[0-9]+[ACGTNacgtn]+|-[0-9]+[ACGTNacgtn]+', '', mpiledup)


        mismatch_dict['A'] = mpileupsnv.count('A')
        mismatch_dict['T'] = mpileupsnv.count('T')
        mismatch_dict['G'] = mpileupsnv.count('G')
        mismatch_dict['C'] = mpileupsnv.count('C')
        trio_depth_counts.update({individual: depths})
        trio_alt_counts.update({individual: mismatch_dict})
        
        

    # now parse the dictionaries into a output string
    father_depth = trio_depth_counts['father']
    mother_depth = trio_depth_counts['mother']
    father_alt_base = [k for k, v in trio_alt_counts['father'].items() if v > 0]
    mother_alt_base = [k for k, v in trio_alt_counts['mother'].items() if v > 0]
    child_alt_base = [k for k, v in trio_alt_counts['child'].items() if v > 0]

#     print(trio_alt_counts)
#     print(father_alt_base, mother_alt_base)
    alt_bases = set(father_alt_base).union(set(mother_alt_base)).union(set(child_alt_base)).intersection(['A', 'C', 'G', 'T'])
#     print(alt_bases)
    snvcount = ''
    for alt_base in alt_bases:
        alt_father = trio_alt_counts['father'][alt_base]
        alt_mother = trio_alt_counts['mother'][alt_base]
        alt_child = trio_alt_counts['child'][alt_base]
        child_depth = trio_depth_counts['child']
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)
        child_vaf = vaf(alt_child, child_depth)
        
        if child_vaf > 0.005:
            snvcount += f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{alt_child}\t{child_depth}\t{child_vaf}\tNA\tNA\t{father_vaf}\t{mother_vaf}\n'

    return snvcount



def get_parent_het_homref_child_count(mpileup_file):
    """parse mpileup results into a table"""
    output_counts_region = mpileup_file +'.counts'
    with open(output_counts_region, 'w') as g:
        g.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\thetero_parent\thomoalt_parent\tfather_vaf\tmother_vaf\n')
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                g.write(parse_mpileup_mtdna(line))

        f.close()
    g.close()

    return output_counts_region
        

def run_mle_rscript(count_table, output_dir, run_mode):
    """run mode can be either optim for quickly running mle, and plot if you want to get mle plot for all possible estimation"""

    cmd = f'{RSCRIPT} {MLE_RSCRIPT} -i {count_table} -o {output_dir} -r {run_mode}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 


def get_paths(path_config):
    """configures the paths to SAMTOOLS AND VARSCAN"""
    with open(path_config) as f:
        path = json.load(f)
    return path['SAMTOOLS'],  path['RSCRIPT'], path['GZIP']


def main():
    global SAMTOOLS, REFERENCE, RSCRIPT, MLE_RSCRIPT, GZIP

    father_bam, mother_bam, child_bam, REFERENCE, mtdna, thread, output_dir, prefix = argument_parser()
    output_dir = os.path.abspath(output_dir)

    # configure paths to executables 
    script_dir = os.path.dirname(os.path.realpath(__file__)) 
    path_config = os.path.join(script_dir, 'path_config.json')

    SAMTOOLS, RSCRIPT, GZIP = get_paths(path_config)

    # path to the MLE Rscript 
    MLE_RSCRIPT = os.path.join(script_dir, 'chimera_likelihood.R')



    # run mpileup parallelly
    tmp_region_dir = os.path.join(output_dir, prefix + '_tmp')
    os.system(f'mkdir -p {tmp_region_dir}')
    arg_list = []
    arg_list.append((father_bam, mother_bam, child_bam, mtdna, tmp_region_dir))

    # run with multithreading
    print('variant calling')
    with mp.Pool(thread) as pool:
        mpileup_files = pool.starmap(mpileup, arg_list)

    # go through mpileup files to parse information
    print('parsing mpileup')
    with mp.Pool(thread) as pool:
        counts_split_files = pool.map(get_parent_het_homref_child_count, mpileup_files)


    print(counts_split_files)


    # combine counts files
    combined_counts = os.path.join(output_dir, f'{prefix}.mtdna.counts')
    count = 0
    with open(combined_counts, 'w') as f:
        for count_file in counts_split_files:
            with open(count_file, 'r') as g:
                if count != 0:
                    next(g) # write header only in the first file, otherwise skip first row
                count += 1
                for line in g:
                    f.write(line)


    # # maximum likelihood estimate
    # print('running MLE')
    # run_mle_rscript(combined_counts, output_dir, run_mode='optim')


if __name__=='__main__':
    main()
