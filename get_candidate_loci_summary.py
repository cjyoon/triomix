

import cyvcf2
import os
import re

import argparse

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--varscan_vcf', required=True, help='')
    parser.add_argument('-f', '--father_id', required=True, help='Readgroup ID of the father in the VCF file')
    parser.add_argument('-m', '--mother_id', required=True, help='Readgroup ID of the mother in the VCF file')
    parser.add_argument('-c', '--child_id', required=True, help='Readgroup ID of the child in the VCF file')
    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('-p', '--prefix', required=True, help='Prefix for the output file')
    parser.add_argument('-d', '--depth_threshold_child', required=True, type=int, default=10, help='threshold of the child at each loci')

    args = vars(parser.parse_args())

    return args['varscan_vcf'], args['output_dir'], args['prefix'], args['father_id'], args['mother_id'], args['child_id'], args['depth_threshold_child']


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



def get_parent_het_homref_child_count(child_id, father_id, mother_id, vcf, prefix, output_dir, depth_threshold_child):
    father_index, mother_index, child_index = identify_trio_index(father_id, mother_id, child_id, vcf)
    count = 0 
    prev_chrom = ''
    chrom_list = get_chromosome_list(vcf)
    autosome_list = [i for i in chrom_list if not re.search(r'Y', i)] # autosome + chrX
    with open(os.path.join(output_dir, f'{prefix}.counts'), 'w') as f:
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

    
            
    print(os.path.join(output_dir, f'{child_id}.vafs'))
    return 0 

def main():
	varscan_vcf, output_dir, prefix, father_id, mother_id, child_id, depth_threshold_child = argument_parser()
	get_parent_het_homref_child_count(child_id, father_id, mother_id, varscan_vcf, prefix, output_dir, depth_threshold_child)
	return 0 

if __name__=='__main__':
    main()


