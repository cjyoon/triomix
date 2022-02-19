import sys
import os
import subprocess
import shlex
import re
import argparse
import multiprocessing as mp
import pysam

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--father', required=True, help='')
    parser.add_argument('-m', '--mother', required=True, help='')
    parser.add_argument('-c', '--child', required=True, help='')
    parser.add_argument('-s', '--sibling', required=False, help='')
    parser.add_argument('-r', '--ratio', nargs="+", required=True, type=float, help='fraction of father, mother, child, sibling ratio in order')

    parser.add_argument('-o', '--output_dir', required=False, default=os.getcwd(), help='Output directory')
    args = vars(parser.parse_args())

    return args['father'], args['mother'], args['child'], args['sibling'], args['output_dir'], args['ratio']

def count_reads(bamfile):
    """count total reads in a bamfile to adjust for the sampling ratio"""
    count= float(0)
    if bamfile != None:
        bam = pysam.AlignmentFile(bamfile)
        for read in bam.fetch():
            count += 1

    return count

def subsample_bam(bamfile, subsample_ratio, subsampled_bam, threads):
    """subsample bam that will be merged later"""
    cmd = f'sambamba view -t {threads} -f bam -h -s {subsample_ratio} {bamfile} > {subsampled_bam}'
    print(cmd)
    os.system(cmd)

    return subsampled_bam

def calculate_adjusted_sampling_ratio(readcount_dict, ratio):
    """calculate thea adjusted sampling ratio that accounts for total read counts of each bam"""

    father_ratio = ratio[0]*readcount_dict['offspring']/readcount_dict['father']
    mother_ratio = ratio[1]*readcount_dict['offspring']/readcount_dict['mother']
    offspring_ratio = ratio[2]
    sibling_ratio = ratio[3]*readcount_dict['offspring']/readcount_dict['sibling']

    return father_ratio, mother_ratio, offspring_ratio, sibling_ratio


def main():
    father, mother, offspring, sibling, output_dir, ratio = argument_parser()
    threads = 4
    os.system(f'mkdir -p {output_dir}')
    if sibling != None:
        with open(os.path.join(output_dir, 'ratio.txt'), 'w') as f:
            f.write(f'case\tfather\tmother\toffspring\tsibling\n{output_dir}\t{ratio[0]}\t{ratio[1]}\t{ratio[2]}\t{ratio[3]}')
    else:
      with open(os.path.join(output_dir, 'ratio.txt'), 'w') as f:
            f.write(f'case\tfather\tmother\toffspring\tsibling\n{output_dir}\t{ratio[0]}\t{ratio[1]}\t{ratio[2]}\t0')

    print(ratio)

    subsampled_father = os.path.join(output_dir, re.sub(r'.bam$', f'_{ratio[0]}.bam', os.path.basename(father)))
    subsampled_mother = os.path.join(output_dir, re.sub(r'.bam$', f'_{ratio[1]}.bam', os.path.basename(mother)))
    subsampled_offspring = os.path.join(output_dir, re.sub(r'.bam$', f'_{ratio[2]}.bam', os.path.basename(offspring)))
    if sibling !=None:
        subsampled_sibling = os.path.join(output_dir, re.sub(r'.bam$', f'_{ratio[3]}.bam', os.path.basename(sibling)))

    merged_bam = os.path.join(output_dir, f'familymix.bam')

    if sibling != None:
        with mp.Pool(4) as pool:
            read_counts = zip(['father', 'mother', 'offspring', 'sibling'], pool.map(count_reads, [father, mother, offspring, sibling]))
    else:
        with mp.Pool(3) as pool:
            read_counts = zip(['father', 'mother', 'offspring'], pool.map(count_reads, [father, mother, offspring]))
       

    readcount_dict= dict()
    for indiv, count in read_counts:
        readcount_dict.update({indiv: count})
 
    with open(os.path.join(output_dir, 'readcounts.txt'), 'w') as f:
        f.write('indiv\ttotal_reads\n')
        for indiv in readcount_dict.keys():
            indiv_readcount = readcount_dict[indiv]
            f.write(f'{indiv}\t{indiv_readcount}\n')

    adjusted_ratio = calculate_adjusted_sampling_ratio(readcount_dict, ratio)
    print(adjusted_ratio)

    with open(os.path.join(output_dir, 'adjusted_ratio.txt'), 'w') as f:
        f.write(f'father\t{adjusted_ratio[0]}\n')
        f.write(f'mother\t{adjusted_ratio[1]}\n')
        f.write(f'offspring\t{adjusted_ratio[2]}\n')
        try:
            f.write(f'sibling\t{adjusted_ratio[3]}\n')
        except:
            pass


    arg_list = []

    arg_list.append((father, adjusted_ratio[0], subsampled_father, threads))
    arg_list.append((mother, adjusted_ratio[1], subsampled_mother, threads))
    arg_list.append((offspring, adjusted_ratio[2], subsampled_offspring, threads))
    if sibling !=None:
        arg_list.append((sibling, adjusted_ratio[3], subsampled_sibling, threads))

    with mp.Pool(4) as pool:
        pool.starmap(subsample_bam, arg_list)

    # merge
    if sibling !=None:
        merge_cmd = f'sambamba merge -t {threads} {merged_bam} {subsampled_father} {subsampled_mother} {subsampled_offspring} {subsampled_sibling}'
    else:
        merge_cmd = f'sambamba merge -t {threads} {merged_bam} {subsampled_father} {subsampled_mother} {subsampled_offspring}'

    print(merge_cmd)

    os.system(merge_cmd)

    print('done')

    
if __name__=='__main__':
    main()
