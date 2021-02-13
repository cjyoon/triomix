"""
chimera scRNAseq decomposition with common SNPs
2020.09.26 cjyoon 
2021.01.07 cjyoon
2021.01.13 cjyoon

"""

import sys
import pysam
import re
import cyvcf2
import multiprocessing as mp
import numpy as np 
import pandas as pd
from collections import Counter, OrderedDict
import argparse
import os


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_rna_bam', required=True, help='Input 10x scRNAseq bam to be deconvoluted')
    parser.add_argument('-hs', '--hetero_snps', nargs='+', required=True, help='heterozygous SNPs in each zygote. Two input files needed if deconvoluting chimera of two zygotes')
    parser.add_argument('-cb', '--cell_barcodes', required=False, help='List of filtered cell barcodes to use')
    parser.add_argument('-m', '--min_bc', required=False, type=int, help='Minimum number of reads per barcode to be used', default=200)
    parser.add_argument('-M', '--max_bc', required=False, type=int, help='Maximum number of reads per barcode to be used', default=90000)
    parser.add_argument('-t', '--thread', required=False, type=int, help='Multithread to utilize', default=2)
    parser.add_argument('-n', '--names', required=True, nargs='+', help='names of each zygote. Will be the in the header of the final output dataframe')

    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-p', '--prefix', required=False, help='output prefix, if not specified, will take the input bam name')
    args = vars(parser.parse_args())

    if args['prefix'] == None:
        args['prefix'] = re.sub(r'.bam$', '', os.path.basename(args['input_rna_bam']))

    return args['input_rna_bam'], args['hetero_snps'], args['cell_barcodes'], args['min_bc'], args['max_bc'], args['output_dir'], args['thread'], args['names'], args['prefix']


def get_reference_contigs(bamfile):
    """get the major chromosomal contigs from the bamfile"""
    bam = pysam.AlignmentFile(bamfile)
    contigs = []
    for contig in bam.references:
        if re.search(r'^chr[0-9XY]+$|^[0-9XY]+$', contig):
            contigs.append(contig)
    bam.close()
    return contigs
           

def get_cb_ub_counts_chrom(rna_bam, chrom):
    """retrieve unique CB, UB barcodes in a 10x RNAseq bam file
    for a specific chromosome
    """
    print(chrom)
    cb_list= []
    ub_list = []
    read_count = 0
    bamfile = pysam.AlignmentFile(rna_bam)
    for read in bamfile.fetch(chrom):
        read_count += 1
        try:
            cb_list.append(read.get_tag('CB'))
            ub_list.append(read.get_tag('UB'))
        except:
            pass
        
    return cb_list, ub_list, read_count

def get_cb_ub_counts(rna_bam, thread=4):
    """get the number of reads for each CB and UB tag in the 10x scRNAseq bam"""
    chromosomes = get_reference_contigs(rna_bam)
    arg_list = zip(np.repeat(rna_bam, len(chromosomes)), chromosomes)
    cb_list = []
    ub_list = []
    read_count = 0
    with mp.Pool(thread) as pool:
        result = pool.starmap(get_cb_ub_counts_chrom, arg_list)
        
    print('summing up')
    for chrom_result in result:
        cb_list += chrom_result[0]
        ub_list += chrom_result[1]
        read_count += chrom_result[2]
        
   
    return Counter(cb_list), Counter(ub_list), read_count

########## these variant detection in read functions come from pztrio_c.py
def read_contains_snv_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """for a given read pair see if this read contains SNV variant or reference"""
    pos = variant_pos - 1
    is_variant_read = 0

    try:
        if (read.get_reference_positions()[0] - pos) * (read.get_reference_positions()[-1] - pos) <= 0:
            position_index_in_read = read.get_reference_positions(
                full_length=True).index(pos)
            base = (read.query_sequence[position_index_in_read]).upper()
            base_quality = read.query_qualities[position_index_in_read]
            if base == variant_alt and base_quality >= 20:
                is_variant_read = True
            elif base == variant_ref and base_quality >= 20:
                is_variant_read = False
            else:
                is_variant_read = None
    #             print(base, base_quality, base_quality >= 20
        else:
            # print('non informative read included somehow')
            raise ValueError
    except:
        is_variant_read = 0  # contains indels, do not use in calculation

    return is_variant_read


def read_contains_insertion_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """for a given read pair see if this read contains Insertion variant or reference"""
    pos = variant_pos - 1
    variant_length = len(variant_alt)
    is_variant_read = 0
    try:
        if (read.get_reference_positions()[0] - pos) * (read.get_reference_positions()[-1] - pos) <= 0:
            insertion_start_index = read.get_reference_positions(
                full_length=True).index(pos) + 1
            insertion_end_index = read.get_reference_positions(
                full_length=True).index(pos) + variant_length
            insertion_base_quality = read.query_qualities[insertion_start_index -
                                                           1:insertion_end_index]
            # only consider variants with insertion base qualities  greater than 20
            if all(i > 10 for i in insertion_base_quality):
                if read.get_reference_positions(full_length=True)[insertion_start_index:insertion_end_index] == [None] * (variant_length-1) and read.query_sequence[insertion_start_index-1:insertion_end_index] == variant_alt:
                    #                 print('read contains variant')
                    is_variant_read = True
                else:
                    is_variant_read = False

        else:
            # print('non informative read included somehow')
            raise ValueError
    except:
        # print('exception')
        is_variant_read = 0

    return is_variant_read


def read_contains_deletion_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """for a given read see if this read contains Deletion variant or reference"""
    pos = variant_pos - 1
    # since deletion variant length is calculated with variant_ref, assumes variant_alt is length 1
    variant_length = len(variant_ref)
    try:
        if (read.get_reference_positions()[0] - pos) * (read.get_reference_positions()[-1] - pos) <= 0:
            deletion_start_index = read.get_reference_positions(
                full_length=True).index(pos)
#             print(read.get_reference_positions(full_length=True))
#             print(deletion_start_index)
#             try:
            if read.get_reference_positions(full_length=True)[deletion_start_index + 1] == variant_pos + variant_length - 1:
                is_variant_read = True
            else:
                is_variant_read = False
#             except:
#                 is_variant_read = 0
        else:
            # print('non informative read included somehow')
            raise ValueError

    except:
        is_variant_read = 0
    return is_variant_read


def read_contains_variant(read, variant_chrom, variant_pos, variant_ref, variant_alt):
    """goes through read covering a position of a variant
    tries to identify those reads that support the provided variant and those that do not
    variant is a cyvcf2 variant class and variant can be SNV, small insertion, or small deletion
    assumes that variants have been decomposed into single alleleic variants
    and also left normalized (ideally using vt)
    2019.02.20 cjyoon
    """
    variant_length = len(variant_alt)
    is_variant_read = False

    # variant is an SNV
    if len(variant_ref) == 1 and len(variant_alt) == 1:
        # print('SNV')
        is_variant_read = read_contains_snv_variant(
            read, variant_chrom, variant_pos, variant_ref, variant_alt)

    # variant is a INSERTION
    elif len(variant_ref) < len(variant_alt):
        # print('INSERTION')
        is_variant_read = read_contains_insertion_variant(
            read, variant_chrom, variant_pos, variant_ref, variant_alt)

    # variant is a DELETION
    # for deletion there is no base quality to check for
    elif len(variant_ref) > len(variant_alt):
        # print('DELETION')
        is_variant_read = read_contains_deletion_variant(
            read, variant_chrom, variant_pos, variant_ref, variant_alt)

    else:
        print(f'{variant_chrom}_{variant_ref}>{variant_alt} is probably not left normalized single allelic record')
        raise ValueError
    return is_variant_read


def parse_variant(variant):
    """parses variant string which is in the following format
    chr1:12345_C>T -> chr1, 12345, C, T
    """
    chrom_pos, ref_alt = variant.split('_')
    chrom, pos = chrom_pos.split(':')
    ref, alt = ref_alt.split('>')
    return chrom, float(pos), ref, alt

    
def parse_input_hetero_snps_file(hetero_snps_file):
    """hetero_snps_file is a dataframe output, which may contain many other columns unnecessary for analysis. 
    This function will take out the column with 'variant' as the column name 
    which contains data in the following format
    chr1:12345_C>T
    """
    df = pd.read_csv(hetero_snps_file, delimiter='\t')

    # now select the variant column
    return sorted(df['variant'])


def filter_cb_min_max(rna_bam, min_bc, max_bc, thread=4):
    """since using all CB in the entire bam is memory inefficient, 
    we will select a subset of CB to use based on the minimum and maximum numbers of reads in the barcode"""
    
    # first get the read counts for each cb and ub. 
    result = get_cb_ub_counts(rna_bam, thread)

    cb_filtered = set()
    for cb, cb_count in result[0].items():
        if cb_count > min_bc and cb_count < max_bc:
            cb_filtered.add(cb)
    return sorted(cb_filtered)

def count_snp_in_scrna(rna_bam_path, cell_barcodes, snp_variants):
    """evaluates the SNP variant status in the RNAseq bam file. 
    Evaluates for each cell barcode whether the SNP variant is present in the read. 
    If the variant is present, then the value in the matrix -> 1
    """
    bamfile = pysam.AlignmentFile(rna_bam_path)
    # initialize the summary dataframe
    snp_barcode = pd.DataFrame(0, index=list(cell_barcodes), columns=snp_variants)
    count = 0 
    for variant in sorted(snp_variants):
        chrom, pos, ref, alt = parse_variant(variant)
        chrom = re.sub(r'^chr', '', chrom)
        if count% 10000 == 0:
            print(variant)
            
        for read in bamfile.fetch(chrom, pos, pos+1):
            if read.has_tag('CB'):
                if (read_contains_variant(read, chrom, pos, ref, alt)) == True:
                    cb = read.get_tag('CB')
                    snp_barcode.loc[cb, variant] = 1
        count += 1

    bamfile.close()
    return snp_barcode

def main():
    scrna_bam, hetero_snps, cell_barcodes, min_bc, max_bc, output_dir, thread, names, prefix = argument_parser()
    with mp.Pool(2) as pool:
        zygote_variants = pool.map(parse_input_hetero_snps_file, hetero_snps)


    # filter for cell barcodes. 
    # Can be a list of barcodes in a file, or can be set as min or max values from argument.
    print('filtering cell barcode based on read counts')
    cb_filtered = filter_cb_min_max(scrna_bam, min_bc, max_bc, thread)


    ###############################################################
    # place holder for cell_barcodes if want to use pre-filtered cell barcode from Seurat
    ###############################################################
    # create an argument list for multiprocessing
    # arg_list = list(zip([scrna_bam]*len(zygote_variants), [cb_filtered]*len(zygote_variants), zygote_variants))
    # print(arg_list)
    # Now for each zygote, create a matrix
    print('Creating matrix for each zygote')

    za_snp_barcode = count_snp_in_scrna(scrna_bam, cb_filtered, zygote_variants[0])
    zb_snp_barcode = count_snp_in_scrna(scrna_bam, cb_filtered, zygote_variants[1])

    # print(snp_barcodes)

    # Now sum up the variant read counts by adding up each row
    print('Summing up variant count for each barcode')
    summed_snp_barcodes = []
    for snp_barcode in [za_snp_barcode, zb_snp_barcode]:
        summed_snp_barcode = snp_barcode.sum(axis=1, skipna=True)
        summed_snp_barcodes.append(summed_snp_barcode)

    # Combine the matrices for the zygotes for the final result
    print('merging results from each zygote')
    summed_snp_barcodes_zygotes = pd.concat(summed_snp_barcodes, axis=1, sort=True)

    # add the column name from the argument
    summed_snp_barcodes_zygotes.columns = names
    #add the cell barcodes as a column
    summed_snp_barcodes_zygotes['cb'] = summed_snp_barcodes_zygotes.index


    # write the final result 
    print('writing results')
    output_file = os.path.join(output_dir, prefix + '.deconv.tsv')
    summed_snp_barcodes_zygotes.to_csv(output_file, index=False, header=True, sep='\t')

    print('done')


    # summary_table = 

if __name__=='__main__':
    main()


