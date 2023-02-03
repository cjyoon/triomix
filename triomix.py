#!/usr/bin/env python3 
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
import random
import pandas as pd
import numpy as np

VERSION='0.0.1'

def argument_parser():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(prog='triomix')
    parser.add_argument('--version', action='version', version='%(prog)s v' + VERSION)
    parser.add_argument('-f', '--father', required=True, help="Father's BAM or CRAM file")
    parser.add_argument('-m', '--mother', required=True, help="Mother's BAM or CRAM file")
    parser.add_argument('-c', '--child', required=True, help="Child's BAM or CRAM file")
    parser.add_argument('-r', '--reference', required=True, help="Reference FASTA file")
    parser.add_argument('-s', '--snp', required=False, default=None, help="Optional list of SNP sites as a BED (or BED.gz) file")
    parser.add_argument('-t', '--thread', required=False, default=1, type=int, help="Multithread to utilize. Default=1")
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory. Default=current working directory')
    parser.add_argument('-p', '--prefix', required=False, default=None, help="prefix for the output file. If not specified, will use the SM tag from the child bam's header")
    parser.add_argument('--runmode', required=False, default='all', choices=['single', 'joint', 'all'], help="Runmode for mle.R script. 'single' assumes only 1 contamination source within family. 'joint' calculates the fraction of all family members jointly. 'all' runs both modes. Default=all")
    parser.add_argument("-u", '--upd', default=1, choices=[0, 1], help="0: mle will filter out vaf=0 or 1 in sites where parental genotypes are homo-ref + homo-alt (GroupA SNPs) 1: mle will identify UPDs which appears as contamination. Default=1")
    parser.add_argument('--parent', required=False, action='store_true', help="Run detection of parental DNA contamination with child's DNA")
    parser.add_argument('-d', '--downsample', required=False, default=0.1, type=float, help="Downsampling for plotting.")

    args = vars(parser.parse_args())

    if args['prefix'] == None:
        args['prefix'] = sampleNameBam(args['child'])


    return args['father'], args['mother'], args['child'], args['reference'], args['snp'], args['thread'], args['output_dir'], args['prefix'], args['runmode'], args['upd'], args['parent'], args['downsample']


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name


def identify_chromosomes(fasta_file):
    """given fasta file, identify the chromosomes that are autosomal"""
    fai_file = fasta_file +  '.fai'
    if(os.path.exists(fai_file)):
        with open(fai_file, 'r') as f:
            for line in f:
                chrom, length, *args = line.split('\t')
                if re.search(r'^chr[0-9XY]+$|^[0-9XY]+$', chrom):
                    yield (chrom, float(length))
    else:
        print(f'There is no index file for {fasta_file}. Exiting...')
        sys.exit(1)

def parse_par_bed_file(par_bed):
    """get start and end position of pseudoautosomal regions from a given bed file"""
    par_exclusion_list = []
    with open(par_bed, 'r') as f:
        for line in f:
            chrom, start, end = line.strip().split('\t')
            par_exclusion_list.append((float(start), float(end)))

    return par_exclusion_list


def split_regions(fasta_file, segment_length):
    """splits chromosome into segment lengths"""
    chr_regions = []
    for chrom, chr_length in identify_chromosomes(fasta_file):
        for i in range(0, max(1, int(chr_length/segment_length))): # fixed here
            start = i*segment_length
            end = (i+1)*segment_length

            if(chr_length-end <= segment_length):
                end = chr_length
            chr_regions.append(f'{chrom}:{start:.0f}-{end:.0f}')
    return chr_regions


def check_gzip_file(file_path):
    """checks if the file is binary or not"""
    cmd = f'file {file_path}'
    execute = subprocess.check_output(shlex.split(cmd)).decode()
    if re.search(r'gzip compressed|gzip compatible', execute):
        return True
    else:
        return False


def check_region_and_snp_bed(region, snp_bed):
    """if a region does not contain anything on the bed varscan output has no variant, it will cause error downstream when merging"""
    region_chromosome, start_end = region.split(':')
    region_start, region_end = start_end.split('-')
    region_start = float(region_start)
    region_end = float(region_end)
    variant_count = 0 

    is_snp_bed_gzip = check_gzip_file(snp_bed)

    # determine gzip compression of snp bed file and use appropriate file handler
    if is_snp_bed_gzip==True:
        f = gzip.open(snp_bed, 'rt')
    else:
        f= open(snp_bed, 'r')


    for line in f:
        chrom, start, end, *args = line.strip().split('\t')
        if chrom == region_chromosome:
            if float(start) > region_start and float(start) < region_end:
                variant_count += 1
            
            if variant_count > 0:
                f.close()
                return True
    # if no variant within the region is found then return False
    f.close()
    return False


def filter_regions_with_snv(region_list, snp_bed):
    """using check_region_and_snp_bed remove regions with no overlaping snp in the region"""
    keep_list = []
    for region in region_list:
        if check_region_and_snp_bed(region, snp_bed):
            keep_list.append(region)
        else:
            print(f'{region} filtered out since it does not have any of the SNPs in the BED file')

    return keep_list

def mpileup(father_bam, mother_bam, child_bam, region, output_dir, snp_bed, prefix):
    """Run mpileup in a defined region of interest"""
    
    region_string = re.sub(r':|-', '_', region)

    # if SNP-BED file is provided 
    if snp_bed != None:
        snp_bed_string = f' -l {snp_bed} '
    else:
        snp_bed_string = ''

    # output_file = os.path.join(output_dir, f'{child_id}_{region_string}.mpileup')
    output_file_compressed = os.path.join(output_dir, f'{prefix}_{region_string}.mpileup.gz')
    # print(output_file_compressed)

    if not os.path.isfile(output_file_compressed):
        cmd = f'{SAMTOOLS} mpileup -B -Q 20 -q 20 {snp_bed_string}-r {region} -f {REFERENCE} {father_bam} {mother_bam} {child_bam} | gzip -f > {output_file_compressed}'
        os.system(cmd)

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


def parse_mpileup(mpileup_line, homoref_sampling_rate):
    """Parse mpileup result into counts string"""
    
    chrom, pos, ref, n_alt_base_father, n_alt_base_mother, n_alt_base_child, \
    father_alt_base, mother_alt_base, child_alt_base, \
    father_depth, mother_depth, child_depth, trio_alt_counts = parse_mpileup_line(mpileup_line)

    if (n_alt_base_father==1 and n_alt_base_mother==0):
        alt_base = ''.join(father_alt_base)
        alt_parent = 'F'
        alt_father = trio_alt_counts['father'][alt_base]
        alt_mother = trio_alt_counts['mother'][alt_base]
        child_alt = trio_alt_counts['child'][alt_base]
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)
    elif (n_alt_base_father==0 and n_alt_base_mother==1):
        alt_base = ''.join(mother_alt_base)
        alt_parent = 'M'
        alt_father = trio_alt_counts['father'][alt_base]
        alt_mother = trio_alt_counts['mother'][alt_base]
        child_alt = trio_alt_counts['child'][alt_base]
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)

    elif (n_alt_base_father==0 and n_alt_base_mother==0):
        alt_base='N' # parents homoref
        alt_parent='NA'
        alt_father = 0
        alt_mother=0
        # if parents are homoref/homoref, any non ref bases are errors
        child_alt = sum(trio_alt_counts['child'].values())
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)

    else:
        father_vaf = -1 # arbitrary to filter out
        mother_vaf = -1 # arbitrary to filter out

        pass
    snvcount = ''
    if (father_vaf > 0.4 and father_vaf < 0.6 and mother_vaf ==0) or (mother_vaf > 0.4 and mother_vaf < 0.6 and father_vaf ==0):
        if father_depth > 10 and mother_depth > 10:
            child_vaf = vaf(child_alt, child_depth)
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\t{alt_parent}\tNA\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'
    elif (father_vaf==0 and mother_vaf ==1) or (father_vaf==1 and mother_vaf==0):
        if father_depth > 10 and mother_depth > 10:
            child_vaf = vaf(child_alt, child_depth)
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\tNA\t{alt_parent}\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'
    elif (father_vaf==0 and mother_vaf==0) and random_sample_selection(homoref_sampling_rate):
        if father_depth > 10 and mother_depth > 10:
            child_vaf = vaf(child_alt, child_depth)
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\tNA\tNA\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'
    return snvcount


def parse_mpileup_line(mpileup_line):
    """due to repeated use, decided to make this into a function"""
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
    child_depth = trio_depth_counts['child']

    father_alt_base = [k for k, v in trio_alt_counts['father'].items() if v > 0]
    mother_alt_base = [k for k, v in trio_alt_counts['mother'].items() if v > 0]
    child_alt_base = [k for k, v in trio_alt_counts['child'].items() if v > 0]

    n_alt_base_father = len(father_alt_base)
    n_alt_base_mother = len(mother_alt_base)
    n_alt_base_child = len(child_alt_base)

    return chrom, pos, ref, n_alt_base_father, n_alt_base_mother, n_alt_base_child, father_alt_base, mother_alt_base, child_alt_base, father_depth, mother_depth, child_depth, trio_alt_counts


def parse_mpileup_child_homoalt(mpileup_line):
    """Parse mpileup result into counts string
    Only find places where child is homoalt"""
    
    chrom, pos, ref, n_alt_base_father, n_alt_base_mother, n_alt_base_child, \
    father_alt_base, mother_alt_base, child_alt_base, \
    father_depth, mother_depth, child_depth, trio_alt_counts = parse_mpileup_line(mpileup_line)

    if (n_alt_base_child==1):
        alt_base = ''.join(child_alt_base)
        alt_father = trio_alt_counts['father'][alt_base]
        alt_mother = trio_alt_counts['mother'][alt_base]
        child_alt = trio_alt_counts['child'][alt_base]
        father_vaf = vaf(alt_father, father_depth)
        mother_vaf = vaf(alt_mother, mother_depth)
        child_vaf = vaf(child_alt, child_depth)
    else:
        father_vaf = -1 # arbitrary to filter out
        mother_vaf = -1 # arbitrary to filter out
        child_vaf = -1 # arbitrary to filter out
        pass
    snvcount = ''
    if child_vaf==1:
        if father_depth > 10 and mother_depth > 10:
            snvcount = f'{chrom}\t{pos:.0f}\t{ref}\t{alt_base}\t{child_alt}\t{child_depth}\t{child_vaf}\t{alt_father}\t{father_depth}\t{father_vaf}\t{alt_mother}\t{mother_depth}\t{mother_vaf}\n'

    return snvcount


def get_counts_childhomoalt(mpileup_file):
    """parse mpileup results into a table"""
    output_counts_region = mpileup_file +'.childhomoalt.counts'
    with open(output_counts_region, 'w') as g:
        g.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\tfather_alt\tfather_depth\tfather_vaf\tmother_alt\tmother_depth\tmother_vaf\n')
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                g.write(parse_mpileup_child_homoalt(line))

        f.close()
    g.close()

    return output_counts_region

def random_sample_selection(sampling_rate):
    """random sampling function, avoid sampling if rate==1"""
    if sampling_rate==1:
        return True
    else:
        prob = random.uniform(0, 1)
        if prob < sampling_rate:
            return True
        else:
            return False


def get_child_count(mpileup_file, homoref_sampling_rate):
    """parse mpileup results into a table"""
    output_counts_region = mpileup_file +'.counts'
    with open(output_counts_region, 'w') as g:
        g.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\thetero_parent\thomoalt_parent\tfather_alt\tfather_depth\tfather_vaf\tmother_alt\tmother_depth\tmother_vaf\n')
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                g.write(parse_mpileup(line, homoref_sampling_rate))

        f.close()
    g.close()

    return output_counts_region
        

def run_mle_rscript(count_table, output_dir, runmode, upd):
    """run mle script"""

    cmd = f'{RSCRIPT} {MLE_RSCRIPT} -i {count_table} -o {output_dir} -r {runmode} -u {upd}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 

def run_mle_parent_rscript(x_combined_counts, output_dir, sexchrom):
    """Run mle script for parent contam by child"""
    cmd = f'{RSCRIPT} {MLE_PARENT_RSCRIPT} -i {x_combined_counts} -o {output_dir} -x {sexchrom}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 

def run_mle_parent_child_homoalt_rscript(child_homoalt_counts, output_dir):
    """Run mle script for parent contam by child"""
    cmd = f'{RSCRIPT} {MLE_PARENT_CHILD_HOMOALT_RSCRIPT} -i {child_homoalt_counts} -o {output_dir}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 


def run_plot_rscript(count_table, output_dir, reference, downsample=0.1):
    """run mle script"""

    cmd = f'{RSCRIPT} {PLOT_RSCRIPT} -i {count_table} -o {output_dir} -r {reference} -d {downsample}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 


def run_plot_parent_rscript(count_table, output_dir, reference, downsample=0.1):
    """run mle script"""

    cmd = f'{RSCRIPT} {PLOT_RSCRIPT_PARENT} -i {count_table} -o {output_dir} -r {reference} -d {downsample}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 


def run_segmentation(count_table, output_dir, segment_length):
    """run segmentation on homo-ref + homo-alt sites for each parental SNPs"""
    cmd = f'{RSCRIPT} {SEGMENTATION_RSCRIPT} -i {count_table} -o {output_dir} -s {segment_length}'
    print(cmd)
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 



def get_paths(path_config):
    """configures the paths to SAMTOOLS AND VARSCAN"""
    with open(path_config) as f:
        path = json.load(f)
    return path['SAMTOOLS'],  path['RSCRIPT'], path['GZIP']


def get_child_count_chrX(mpileup_file, sexchrom, par_regions):
    """parse mpileup results into a table"""
    output_counts_region = mpileup_file +'.chrX.counts'
    with open(output_counts_region, 'w') as g:
        g.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\tchrX_group\tfather_alt\tfather_depth\tfather_vaf\tmother_alt\tmother_depth\tmother_vaf\n')
        with gzip.open(mpileup_file, 'rt') as f:
            for line in f:
                g.write(parse_mpileup_chrX(line, sexchrom, par_regions))

        f.close()
    g.close()

    return output_counts_region

def position_pseudoautosomal(pos, par_list):
    """True if position is pseudoautosomal, False if position is not pseudoautosomal"""
    status = False

    # if position is inside any of the par list, then it is in PAR region
    for start, end in par_list:
        if (pos-start) * (pos-end) < 0:
            status = True


    return status



def x_to_autosome_ratio(counts_file, individual):
    """measures the depth of autosome to chrX to determine the sex of each individual
    individual can be 'father', 'mother', or '' (child)
    """
    df = pd.read_csv(counts_file, sep='\t')
    if individual!='child':
        individual = individual + '_'
    else:
        individual = ''

    x_depth = np.nanmean(df[df['chrom'].str.contains('^X$|^chrX$')==True][[f'{individual}depth']])
    autosome_depth = np.nanmean(df[df['chrom'].str.contains('^X$|^chrX$')==False][[f'{individual}depth']])
    x_to_auto = float(x_depth/autosome_depth)

    return x_to_auto

def sexchrom_ratio(counts_file):
    
    child_x_ratio = x_to_autosome_ratio(counts_file, 'child')
    father_x_ratio = x_to_autosome_ratio(counts_file, 'father')
    mother_x_ratio = x_to_autosome_ratio(counts_file, 'mother')

    return child_x_ratio, father_x_ratio, mother_x_ratio


def combine_count_files(file_list, output_file):
    """combines multiple split up count file into one single file"""
    count = 0
    with open(output_file, 'w') as f:
        for afile in file_list:
            with open(afile, 'r') as g:
                if count != 0:
                    next(g) # write header only in the first file, otherwise skip first row
                count += 1
                for line in g:
                    f.write(line)
    return output_file

def main():
    global SAMTOOLS, REFERENCE, RSCRIPT, MLE_RSCRIPT, GZIP, PLOT_RSCRIPT, SEGMENTATION_RSCRIPT, MLE_PARENT_RSCRIPT, MLE_PARENT_CHILD_HOMOALT_RSCRIPT, PLOT_RSCRIPT_PARENT, VBID


    father_bam, mother_bam, child_bam, REFERENCE, snp_bed, thread, output_dir, prefix, runmode, upd, parent, downsample, *args  = argument_parser()
    output_dir = os.path.abspath(output_dir)

    # configure paths to executables 
    script_dir = os.path.dirname(os.path.realpath(__file__)) 
    path_config = os.path.join(script_dir, 'path_config.json')

    SAMTOOLS, RSCRIPT, GZIP = get_paths(path_config)

    # path to the MLE Rscript 
    MLE_RSCRIPT = os.path.join(script_dir, 'mle.R')

    # path to the MLE Rscript for parent contamination by child
    MLE_PARENT_RSCRIPT = os.path.join(script_dir, 'mle_parent.R')
    MLE_PARENT_CHILD_HOMOALT_RSCRIPT = os.path.join(script_dir, 'mle_parent.R')

    # path to the plot Rscript 
    PLOT_RSCRIPT = os.path.join(script_dir, 'plot_variant.R')

    # path to the plot parent Rscript 
    PLOT_RSCRIPT_PARENT = os.path.join(script_dir, 'plot_variant_parent.R')

    # path to segmentation Rscript
    SEGMENTATION_RSCRIPT = os.path.join(script_dir, 'upd_segmentation.R')


    # split up regions
    segment_length = 50000000
    region_splits = split_regions(REFERENCE, segment_length)

    # if SNP bed is supplied, filter out those regions where no variant is found. 
    if snp_bed != None:
        region_splits= filter_regions_with_snv(region_splits, snp_bed)
    else:
        pass

    # run mpileup parallelly
    tmp_region_dir = os.path.join(output_dir, prefix + '_tmp')
    os.system(f'mkdir -p {tmp_region_dir}')
    arg_list = []
    for region in region_splits:
        arg_list.append((father_bam, mother_bam, child_bam, region, tmp_region_dir, snp_bed, prefix))

    # run with multithreading
    print('variant calling')
    with mp.Pool(thread) as pool:
        mpileup_files = pool.starmap(mpileup, arg_list)

    # go through mpileup files to parse information
    print('parsing mpileup')
    # prepare arglist
    arg_list = []
    if not snp_bed is None:
        homoref_sampling_rate = 1
    else:
        homoref_sampling_rate = 0.01

    for mpileup_f in mpileup_files:
        arg_list.append((mpileup_f, homoref_sampling_rate))


    with mp.Pool(thread) as pool:
        counts_split_files = pool.starmap(get_child_count, arg_list)


    # print(counts_split_files)


    # combine counts files

    combined_counts = os.path.join(output_dir, f'{prefix}.child.counts')


    # combine the count files for each split mpileup
    combine_count_files(counts_split_files, combined_counts)

    # maximum likelihood estimate
    print('running MLE')
    run_mle_rscript(combined_counts, output_dir, runmode, upd)

    # plot variants
    print('plotting variants')
    run_plot_rscript(combined_counts, output_dir, REFERENCE, downsample=downsample)

    # upd segmentation
    print('upd segmentation')
    upd_segment_length = 1000000 # 1mb
    run_segmentation(combined_counts, output_dir, segment_length=upd_segment_length)


    # get the chrX to autosome depth ratio. 
    child_x_ratio, father_x_ratio, mother_x_ratio = sexchrom_ratio(combined_counts)
    with open(os.path.join(output_dir, prefix + '.x2a.depth.tsv'), 'w') as f:
        f.write(f'child\t{child_x_ratio:.3f}\nfather\t{father_x_ratio:.3f}\nmother\t{mother_x_ratio:.3f}')


    # run parent DNA contamination if --parent is set
    print(f'Parent mode {parent}')
    if parent == True: 
        # AUTOSOME
        with mp.Pool(thread) as pool:
            parent_counts_split_files = pool.map(get_counts_childhomoalt, mpileup_files)

        parent_counts = os.path.join(output_dir, f'{prefix}.parent.counts')
        combine_count_files(parent_counts_split_files, parent_counts)
        print('running MLE on parent contamination')
        run_mle_parent_child_homoalt_rscript(parent_counts, output_dir)
        print('plotting parent variants')
        run_plot_parent_rscript(parent_counts, output_dir, REFERENCE, downsample=downsample)




    print('done')


if __name__=='__main__':
    python_version = sys.version_info[0] + 0.1*sys.version_info[1]
    if python_version < 3.5:
        print('Triomix requires python version 3.5 or later.')
        print('You can specify the python interpreter such as python triomix -f father.bam -m mother.bam -c child.bam')
        print('Exiting...')
    main()
