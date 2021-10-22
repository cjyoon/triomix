import sys
import os
import subprocess
import shlex
import re
import argparse
import json
import multiprocessing as mp
import pysam
import cyvcf2
import gzip

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--father', required=True, help="Father's bam file")
    parser.add_argument('-m', '--mother', required=True, help="Mother's bam file")
    parser.add_argument('-c', '--child', required=True, help="Child's bam file")
    parser.add_argument('-r', '--reference', required=True, help="Reference fasta file")
    parser.add_argument('-s', '--snp', required=False, default=None, help="Optional list of SNP sites as a BED file")
    parser.add_argument('-t', '--thread', required=False, default=1, type=int, help="Multithread to utilize")
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-p', '--prefix', required=False, default=None, help="prefix for the output file. If not specified, will use the SM tag from the child's bam")
    args = vars(parser.parse_args())

    if args['prefix'] == None:
        args['prefix'] = sampleNameBam(args['child'])
    return args['father'], args['mother'], args['child'], args['reference'], args['snp'], args['thread'], args['output_dir'], args['prefix']


def sampleNameBam(bamFile):
    """get @RG SM: information as sample name from BAM header"""
    bam = pysam.AlignmentFile(bamFile)
    name = bam.header['RG'][0]['SM']
    return name


def identify_autosomal_chromosomes(fasta_file):
    """given fasta file, identify the chromosomes that are autosomal"""
    fai_file = fasta_file +  '.fai'
    if(os.path.exists(fai_file)):
        with open(fai_file, 'r') as f:
            for line in f:
                chrom, length, *args = line.split('\t')
                if re.search(r'^chr[0-9]+$|^[0-9]+$', chrom):
                    yield (chrom, float(length))
    else:
        print(f'There is no index file for {fasta_file}. Exiting...')
        sys.exit(1)


def split_regions(fasta_file, segment_length):
    """splits chromosome into segment lengths"""
    chr_regions = []
    for chrom, chr_length in identify_autosomal_chromosomes(fasta_file):
        for i in range(0, int(chr_length/segment_length)):
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
    if re.search(r'gzip compressed', execute):
        return True
    else:
        return False


def write_sample_list(output_dir, father_bam, mother_bam, child_bam, prefix):
    """write the sample name list to be used for varscan vcf labeling
    order always father, mother, child"""
    os.system(f'mkdir -p {output_dir}/samplelist')
    output_file = os.path.join(output_dir, f'samplelist/{prefix}.samplelist')
    with open(output_file, 'w') as f:
        for individual in [father_bam, mother_bam, child_bam]:
            individual_id = sampleNameBam(individual)
            f.write(individual_id + '\n')
    return output_file


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


def varscan_mpileup2snv(father_bam, mother_bam, child_bam, sample_list_file, region, output_dir, snp_bed):
    """Run varscan mpileup2snv in a defined region of interest"""
    
    child_id = sampleNameBam(child_bam)
    region_string = re.sub(r':|-', '_', region)
    # print(region_string)

    # if SNP-BED file is provided 
    if snp_bed != None:
        snp_bed_string = f' -l {snp_bed} '
    else:
        snp_bed_string = ''

    output_file = os.path.join(output_dir, f'{child_id}_{region_string}.varscan.snv.vcf')
    output_file_compressed = os.path.join(output_dir, f'{child_id}_{region_string}.varscan.snv.vcf.gz')
    f = open(output_file, 'w')

    mpileup_cmd = f'{SAMTOOLS} mpileup -B -Q 20 -q 20 {snp_bed_string}-r {region} -f {REFERENCE} {father_bam} {mother_bam} {child_bam}'
    mpileup_execute = subprocess.Popen(shlex.split(mpileup_cmd), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    varscan_cmd = f'{JAVA} -jar {VARSCAN} mpileup2snp --min-coverage 10 --mean-reads2 2  --min-var-freq 0.01 --p-value 0.99 --output-vcf --strand-filter 1 --vcf-sample-list {sample_list_file}' 
    varscan_execute = subprocess.Popen(shlex.split(varscan_cmd), stdin=mpileup_execute.stdout, stdout=f, stderr=subprocess.DEVNULL)
    
    mpileup_execute.stdout.close()
    output = varscan_execute.communicate()[0]
    output_file_compressed = bgzip_tabix(output_file)

    return  output_file_compressed


def combine_vcf_files(split_vcfs, output_dir, child_bam):
    """combine vcf files from parallel jobs into one vcf file"""
    child_id = sampleNameBam(child_bam)
    output_file = os.path.join(output_dir, f'{child_id}.varscan.snv.vcf.gz')
    split_vcfs_string =  ' '.join(split_vcfs)
    cmd = f'{BCFTOOLS} concat {split_vcfs_string} -o {output_file} -O z -a'

    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return output_file


def bgzip_tabix(vcf_file):
    """compress and tabix vcf file"""
    output_file_compressed = re.sub(r'.vcf$', '.vcf.gz', vcf_file)

    bgzip_cmd =  f'{BGZIP} -f {vcf_file}'
    bgzip_execute = subprocess.Popen(shlex.split(bgzip_cmd))
    bgzip_execute.wait()

    tabix_cmd = f'{TABIX} -p vcf {output_file_compressed}' 
    tabix_execution = subprocess.Popen(shlex.split(tabix_cmd))
    tabix_execution.wait()

    return output_file_compressed


def vt_decompose_normalize(vcf, output_dir):
    """decompose and normalize using VT"""
    output_vcf = os.path.join(output_dir, re.sub(r'.vcf.gz$', '.vtdcn.vcf', os.path.basename(vcf)))
    f = open(output_vcf, 'w')
    decompose_cmd = f'{VT} decompose -s  {vcf}'
    normalize_cmd = f'{VT} normalize -n -r {REFERENCE} - '
    decompose_execute = subprocess.Popen(shlex.split(decompose_cmd), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    normalize_execute = subprocess.Popen(shlex.split(normalize_cmd), stdin=decompose_execute.stdout, stdout=f, stderr=subprocess.DEVNULL)
    decompose_execute.stdout.close()
    output = normalize_execute.communicate()[0]

    f.close()
    output_vcf_compressed = bgzip_tabix(output_vcf)
    return output_vcf_compressed

def remove_non_acgt(vcf, output_dir):
    """remove if ref contains a non acgt bases"""
    output_vcf = os.path.join(output_dir, re.sub(r'.vcf.gz$', '.acgt.vcf', os.path.basename(vcf)))
    vcf_handle = cyvcf2.VCF(vcf)
    output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
    for variant in vcf_handle:
        if re.search(r'[ATGC]+', variant.REF) and re.search(r'[ATGC]+', variant.ALT[0]):
            output_handle.write_record(variant)

    vcf_handle.close()
    output_handle.close()

    return bgzip_tabix(output_vcf)


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

def identify_trio_index(father_id, mother_id, child_id, vcf):
    """get the index positions of father mother child trio in the 
    given vcf file. 
    Infer this from the VCF header

    """
    vcf_handle = cyvcf2.VCF(vcf)

    try:
        # start indexing from 9 since first 9 is required vcf columns
        father_index = [i for i, j in enumerate(
            vcf_handle.samples) if j == father_id][0]
        mother_index = [i for i, j in enumerate(
            vcf_handle.samples) if j == mother_id][0]
        child_index = [i for i, j in enumerate(
            vcf_handle.samples) if j == child_id][0]
        return father_index, mother_index, child_index

    except IndexError:
        print(f'Given {vcf} file does not contain the family members {father_id}, {mother_id}, {child_id}')
        sys.exit(1)


def natural_sort(l):
    def convert(text): return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): return [convert(c)
                                   for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def get_chromosome_list(vcf):
    """get all chromosome list that appears in a vcf that is not an alternative contig"""
    chromosomes = set()
    vcf_handle = cyvcf2.VCF(vcf)
    for i in vcf_handle.header_iter():
        if i['HeaderType'] == 'CONTIG':
            chromosomes.add(i['ID'])
        
    chromosomes_list = [chrom for chrom in chromosomes if re.search(r'^chr[0-9XY]+$|^[0-9XY]+$', chrom)]

    return natural_sort(list(chromosomes_list))


def get_parent_het_homref_child_count(child_bam, father_bam, mother_bam , vcf, prefix, output_dir, depth_threshold_child):

    father_id = sampleNameBam(father_bam)
    mother_id = sampleNameBam(mother_bam)
    child_id = sampleNameBam(child_bam)

    father_index, mother_index, child_index = identify_trio_index(father_id, mother_id, child_id, vcf)
    count = 0 
    prev_chrom = ''
    chrom_list = get_chromosome_list(vcf)
    autosome_list = [i for i in chrom_list if not re.search(r'Y', i)] # autosome + chrX
    count_table_path = os.path.join(output_dir, f'{prefix}.counts')
    with open(count_table_path, 'w') as f:
        f.write('chrom\tpos\trefbase\taltbase\talt\tdepth\tvaf\thetero_parent\n')

        for variant in cyvcf2.VCF(vcf):
            if  variant.CHROM in autosome_list:
                current_chrom = variant.CHROM
                father_alt = count_int(variant.format('AD')[father_index][0])
                mother_alt = count_int(variant.format('AD')[mother_index][0])
                child_alt = count_int(variant.format('AD')[child_index][0])

                father_ref = count_int(variant.format('RD')[father_index][0])
                mother_ref = count_int(variant.format('RD')[mother_index][0])
                child_ref = count_int(variant.format('RD')[child_index][0])

                father_depth = int(father_ref) + int(father_alt)
                mother_depth = int(mother_ref) + int(mother_alt)
                child_depth = int(child_ref) + int(child_alt)

                father_vaf = vaf(father_alt, father_depth)
                mother_vaf = vaf(mother_alt, mother_depth)
                child_vaf = vaf(child_alt, child_depth)
                # if prev_chrom != current_chrom:
                #     print(current_chrom)
                prev_chrom = current_chrom
                if father_depth > 20 and mother_depth > 20 and child_depth > depth_threshold_child:
                    if (father_vaf > 0.4 and father_vaf < 0.6 and mother_vaf< 0.01):
                        f.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{child_alt}\t{child_depth}\t{child_vaf}\tF\n')

                    elif (mother_vaf > 0.4 and mother_vaf < 0.6 and father_vaf < 0.01):
                        f.write(f'{variant.CHROM}\t{variant.POS}\t{variant.REF}\t{variant.ALT[0]}\t{child_alt}\t{child_depth}\t{child_vaf}\tM\n')
                    else:
                        pass
    return count_table_path

def run_mle_rscript(count_table, output_dir, run_mode):
    """run mode can be either optim for quickly running mle, and plot if you want to get mle plot for all possible estimation"""

    cmd = f'{RSCRIPT} {MLE_RSCRIPT} -i {count_table} -o {output_dir} -r {run_mode}'
    execute = subprocess.Popen(shlex.split(cmd))
    execute.wait()

    return 0 


def get_paths(path_config):
    """configures the paths to SAMTOOLS AND VARSCAN"""
    with open(path_config) as f:
        path = json.load(f)
    return path['SAMTOOLS'], path['JAVA'], os.path.join(os.path.dirname(os.path.realpath(__file__)), path['VARSCAN']), path['BGZIP'], path['TABIX'], path['BCFTOOLS'], path['VT'], path['RSCRIPT']


def main():
    global SAMTOOLS, JAVA, VARSCAN, BGZIP, TABIX, BCFTOOLS, VT, REFERENCE, RSCRIPT, MLE_RSCRIPT

    father_bam, mother_bam, child_bam, REFERENCE, snp_bed, thread, output_dir, prefix = argument_parser()
    output_dir = os.path.abspath(output_dir)

    # write sample name for varscan
    sample_list_file = write_sample_list(output_dir, father_bam, mother_bam, child_bam, prefix)

    # configure paths to executables 
    script_dir = os.path.dirname(os.path.realpath(__file__)) 
    path_config = os.path.join(script_dir, 'path_config.json')

    SAMTOOLS, JAVA, VARSCAN, BGZIP, TABIX, BCFTOOLS, VT, RSCRIPT = get_paths(path_config)

    # path to the MLE Rscript 
    MLE_RSCRIPT = os.path.join(script_dir, 'chimera_likelihood.R')

    # split up regions
    segment_length = 50000000
    region_splits = split_regions(REFERENCE, segment_length)

    # if SNP bed is supplied, filter out those regions where no variant is found. 
    if snp_bed != None:
        region_splits= filter_regions_with_snv(region_splits, snp_bed)
    else:
        pass

    # run varscan parallelly
    tmp_region_dir = os.path.join(output_dir, prefix + '_tmp')
    os.system(f'mkdir -p {tmp_region_dir}')
    arg_list = []
    for region in region_splits:
        arg_list.append((father_bam, mother_bam, child_bam, sample_list_file, region, tmp_region_dir, snp_bed))

    # run with multithreading
    print('variant calling')
    with mp.Pool(thread) as pool:
        split_vcfs = pool.starmap(varscan_mpileup2snv, arg_list)

    # combine varscan result
    print('combining variant calls')
    combined_vcf = combine_vcf_files(split_vcfs, output_dir, child_bam)
    
    # VT normalize variants
    print('normalizing variant calls with vt')
    vt_vcf = vt_decompose_normalize(combined_vcf, output_dir)

    # remove non acgt results
    print('removing loci with non ACGT bases')
    acgt_vcf = remove_non_acgt(vt_vcf, output_dir)

    # create count summary table
    print('creating count table')
    count_table = get_parent_het_homref_child_count(child_bam, father_bam, mother_bam, acgt_vcf, prefix, output_dir, depth_threshold_child=10)

    # maximum likelihood estimate
    print('running MLE')
    run_mle_rscript(count_table, output_dir, run_mode='optim')


if __name__=='__main__':
    main()
