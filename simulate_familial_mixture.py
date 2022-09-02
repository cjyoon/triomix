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
    """count total autosomal reads in a bamfile to adjust for the sampling ratio"""
    count= float(0)
    if bamfile != None:
        bam = pysam.AlignmentFile(bamfile)
        autosome_list = []
        for contig in (bam.header['SQ']):
            if re.search(r'^chr[0-9]+$|^[0-9]+$', contig['SN']):
                autosome_list.append(contig['SN'])
        print(autosome_list)
        for autosome in autosome_list:
            print(autosome)
            for read in bam.fetch(autosome):
                count += 1

    return count



def subsample_bam(bamfile, subsample_ratio, subsampled_bam, threads):
    """subsample bam that will be merged later"""
    cmd = f'samtools view -@ {threads-1} --subsample {subsample_ratio} -hbo {subsampled_bam} {bamfile}'
    print(cmd)
    os.system(cmd)

    return subsampled_bam

def calculate_adjusted_sampling_ratio(readcount_dict, ratio):
    """calculate thea adjusted sampling ratio that accounts for total read counts of each bam"""
    # use the sample with minimal # of reads as a base
    min_read_indiv = min(readcount_dict, key=readcount_dict.get)
    print(min_read_indiv)
    min_read = readcount_dict[min_read_indiv]
    
    father_ratio = ratio[0]*min_read/readcount_dict['father']
    mother_ratio = ratio[1]*min_read/readcount_dict['mother']
    offspring_ratio = ratio[2]*min_read/readcount_dict['offspring']
    sibling_ratio = ratio[3]*min_read/readcount_dict['sibling']
        
    return father_ratio, mother_ratio, offspring_ratio, sibling_ratio


def main():
    father, mother, offspring, sibling, output_dir, ratio = argument_parser()
    threads = 4
    merged_bam = os.path.join(output_dir, f'familymix.bam')
    if not os.path.exists(merged_bam):
        os.system(f'mkdir -p {output_dir}')
        if sibling != None:
            with open(os.path.join(output_dir, 'ratio.txt'), 'w') as f:
                f.write(f'case\tfather\tmother\toffspring\tsibling\n{output_dir}\t{ratio[0]}\t{ratio[1]}\t{ratio[2]}\t{ratio[3]}')
        else:
          with open(os.path.join(output_dir, 'ratio.txt'), 'w') as f:
                f.write(f'case\tfather\tmother\toffspring\tsibling\n{output_dir}\t{ratio[0]}\t{ratio[1]}\t{ratio[2]}\t0')

        print(ratio)

        subsampled_father = os.path.join(output_dir, re.sub(r'.bam$|.cram$', f'_{ratio[0]}.bam', os.path.basename(father)))
        subsampled_mother = os.path.join(output_dir, re.sub(r'.bam$|.cram$', f'_{ratio[1]}.bam', os.path.basename(mother)))
        subsampled_offspring = os.path.join(output_dir, re.sub(r'.bam$|.cram$', f'_{ratio[2]}.bam', os.path.basename(offspring)))
        if sibling !=None:
            subsampled_sibling = os.path.join(output_dir, re.sub(r'.bam$|.cram$', f'_{ratio[3]}.bam', os.path.basename(sibling)))


        if sibling != None:
            with mp.Pool(4) as pool:
                family_members = ['father', 'mother', 'offspring', 'sibling']
                family_bam = [father, mother, offspring, sibling]
                read_counts = zip(family_members, pool.map(count_reads, [father, mother, offspring, sibling]))
                output_files = [subsampled_father, subsampled_mother, subsampled_offspring, subsampled_sibling]

        else:
            with mp.Pool(3) as pool:
                family_members = ['father', 'mother', 'offspring']
                family_bam = [father, mother, offspring]
                read_counts = zip(family_members, pool.map(count_reads, [father, mother, offspring]))
                output_files = [subsampled_father, subsampled_mother, subsampled_offspring]


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


        for indiv, adj_ratio, subsampled_bam in zip(family_bam, adjusted_ratio, output_files):
            if adj_ratio > 0: # samtools view -s cannot subsample 0 reads, unlike sambamba
                arg_list.append((indiv, adj_ratio, subsampled_bam, threads))


        with mp.Pool(4) as pool:
            subsampled_bam_list = pool.starmap(subsample_bam, arg_list)

        print(subsampled_bam_list)
        subsampled_bam_list_string = ' '.join(subsampled_bam_list)
        # merge

        merge_cmd = f'samtools merge -@ {threads-1} -o {merged_bam} {subsampled_bam_list_string}'

        print(merge_cmd)

        os.system(merge_cmd)

   
    cmd = f'samtools index -@ {threads-1} {merged_bam}'
    os.system(cmd)


    print('done')

    
if __name__=='__main__':
    main()
