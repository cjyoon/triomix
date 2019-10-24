import cyvcf2
import re
import sys
import subprocess
import shlex

def bgzip_tabix(vcf_file):
    '''compress and index vcf file'''
    indexed_vcf = vcf_file + '.gz'
    cmd = f'bgzip -f {vcf_file}'
    compress = subprocess.Popen(shlex.split(cmd))
    compress.wait()
    cmd = f'tabix -f -p vcf {indexed_vcf}'
    tabix = subprocess.Popen(shlex.split(cmd))
    tabix.wait()
    return indexed_vcf


input_vcf = sys.argv[1]        
output_vcf = re.sub(r'.vcf.gz$|vcf$', '.acgt.vcf', input_vcf)
vcf_handle = cyvcf2.VCF(input_vcf)
output_handle = cyvcf2.Writer(output_vcf, vcf_handle)
for variant in vcf_handle:
    if re.search(r'[ATGC]+', variant.REF) and re.search(r'[ATGC]+', variant.ALT[0]):
        output_handle.write_record(variant)
        
vcf_handle.close()
output_handle.close()

bgzip_tabix(output_vcf )
