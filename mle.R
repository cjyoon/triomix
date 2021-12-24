# infer the most likely chimeric ratio from data
# 2019.08.07 cjyoon
# 2021.12.19 cjyoon added/fixed parent mix
# 2021.12.23 cjyoon joint analysis of offspring, sibling, mother, and father.
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(bbmle))
options(warn=-1)


option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input VAF sites at HET/HOMREF parental sites", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


###############################################################################
##### vectorized likelihood calculations and negative-log-likelihood function 

###############################################################################
# simple mix of sibling only
likelihood_sibling<-function(alt, depth, ratio){
  p_00 = dbinom(alt, depth, prob=0)
  p_10 = dbinom(alt, depth, prob=ratio/2)
  p_01 = dbinom(alt, depth, prob= 0.5*(1-ratio))
  p_11 = dbinom(alt, depth, prob=0.5)
  prob_total = (p_00 + p_10 + p_01 + p_11 + 1e-100)/4
  return(log10(prob_total))
}

sibling_nll<-function(df, par){
  if(par>0){
    return(-sum(likelihood_sibling(df$alt, df$depth, par)))
  }else{
    return(-10000000000)
  }
}

###############################################################################
# parent mix only
likelihood_parent_homohomo<-function(alt, depth, homo_parent, parent_diff){
  if(homo_parent=='F'){
    altfraction=min(1, max(0, (1-parent_diff)/2))
    p = dbinom(alt, depth, prob=altfraction)
    prob_total = (p + 1e-100)
    return(log10(prob_total))
  }else{
    altfraction=min(1, max(0, (1+parent_diff)/2))
    p = dbinom(alt, depth, prob=altfraction)
    prob_total = (p + 1e-100)
    return(log10(prob_total))
  }
}
likelihood_parent_homohomo_vec = Vectorize(likelihood_parent_homohomo) #vectorized due to 'ifelse'

parent_diff_nll <-function(df, par){
  if(par[1] >= -1 | par[1] <=1){
    parent_nll_result = -sum(likelihood_parent_homohomo_vec(df$alt, df$depth, df$homoalt_parent, par))
    return(parent_nll_result)
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}

###############################################################################
# offspring, sibling, father, mother mix
likelihood_family_mix_father_altsnp<-function(alt, depth, k, x, z){
      # k is a constant (x-y = k)
      # proaband het, sibling homoref
      alt_fraction_p10 = min(1, max(0, (1-x-z)/2))
      p10 = dbinom(alt, depth, prob=alt_fraction_p10)
      # proband het, sibling het
      alt_fraction_p11 = min(1, max(0, (1-x)/2))
      p11 = dbinom(alt, depth, prob=alt_fraction_p11)
      
      # proband homoref sibling homoref
      alt_fraction_p00 = min(1, max(0, (x-k)/2))
      p00 = dbinom(alt, depth, prob=alt_fraction_p00)
      
      # proband homoref, sibling het
      alt_fraction_p01 = min(1, max(0, (x-k+z)/2))
      p01 = dbinom(alt, depth, prob=alt_fraction_p01)

    p = p00 + p01 + p10 + p11
  
  prob_total = (p + 1e-100)
  return(log10(prob_total))
}

likelihood_family_mix_mother_altsnp<-function(alt, depth, k, x, z){
  # k is a constant (x-y = k)
  # proaband het, sibling homoref
  alt_fraction_p10 = min(1, max(0, (1-x+k-z)/2))
  p10 = dbinom(alt, depth, prob=alt_fraction_p10)
  # proband het, sibling het
  alt_fraction_p11 = min(1, max(0, (1-x+k)/2))
  p11 = dbinom(alt, depth, prob=alt_fraction_p11)
  
  # proband homoref sibling homoref
  alt_fraction_p00 = min(1, max(0, (x)/2))
  p00 = dbinom(alt, depth, prob=alt_fraction_p00)
  
  # proband homoref, sibling het
  alt_fraction_p01 = min(1, max(0, (x+z)/2))
  p01 = dbinom(alt, depth, prob=alt_fraction_p01)
  
  p = p00 + p01 + p10 + p11
  
  prob_total = (p + 1e-100)
  return(log10(prob_total))
}

family_mix_nll<-function(df, k, x, z){
  if((x >= k & z >= 0) & (2*x+z <=1 + k)) {
    father_het_df = df %>% filter(hetero_parent=='F')
    mother_het_df = df %>% filter(hetero_parent=='M')
    
    father_nll_result = -sum(likelihood_family_mix_father_altsnp(father_het_df$alt, father_het_df$depth, k, x, z))
    mother_nll_result = -sum(likelihood_family_mix_mother_altsnp(mother_het_df$alt, mother_het_df$depth, k, x, z))
    nll_sum = father_nll_result + mother_nll_result
    return(nll_sum)
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}



######################################
# Read in input file and specify output files

if(opt$output_dir=='.' || is.na(opt$output_dir)){
  opt$output_dir = getwd() 
}
# create output directory if it doesn't exist
if(!dir.exists(opt$output_dir)){
  system(paste0('mkdir -p ', opt$output_dir))
}
# argument parsing
input_file = opt$input_file
output_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.mle.pdf')))
summary_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.summary.tsv')))


# read in count table
df <- read_delim(input_file, delim='\t', na=c('', 'NA'), col_types=cols(chrom=col_character(), pos=col_double(), alt=col_double(), depth=col_double(), vaf=col_double(), father_vaf=col_double(), mother_vaf=col_double()))

# filter sex chromosomes
df <- df %>% filter(!str_detect(chrom, pattern = 'chrX|X|chrY|Y')) # remove chrX for calculating MLE estimates

# filter out those regions where both parents are homoref/homoref for Mendelian error rate calculation
df_homorefhomoref <- df %>% filter(is.na(hetero_parent) & is.na(homoalt_parent))
mendelian_error_rate = sum(df_homorefhomoref$alt)/ sum(df_homorefhomoref$depth) * 2 # multiplied by 2 since diploid genome

# identify homo alt homo ref
df_homoalthomoref <- df %>% filter(is.na(hetero_parent) & !is.na(homoalt_parent))


# filter out those regions where both parents are homoref/homoref
df_homoref_het <- df %>% filter(!is.na(hetero_parent) & is.na(homoalt_parent))


# Estimate MLE chimera ratio and confidence interval using mle function
# calcaulate using homo ref + homo_alt loci
df_homoalthomoref <- df %>% filter(is.na(hetero_parent) & !is.na(homoalt_parent))
mother_homoalt = df_homoalthomoref %>% filter(homoalt_parent=='M')
father_homoalt = df_homoalthomoref %>% filter(homoalt_parent=='F')

parent_diff_result = mle2(parent_diff_nll, start=list(par=0), lower=list(par=-1), upper=list(par=1), data=list(df=df_homoalthomoref), method='L-BFGS-B')
parent_diff_value = coef(parent_diff_result)[['par']]

# using the differnce between (mother-father=parent_diff) use this as a constraint, then find values for mother(x) and sibling(z) 
mother_sibling_mle = mle2(family_mix_nll, start=list(x=max(parent_diff_value, 0), z=0), lower=list(x=max(parent_diff_value, 0), z=0), upper=list(x=(1+parent_diff_value)/2, z=0.5), fixed=list(k=parent_diff_value), data=list(df=df_homoref_het), method='L-BFGS-B')

mother_fraction = coef(mother_sibling_mle)[['x']]
sibling_fraction = coef(mother_sibling_mle)[['z']]
father_fraction = mother_fraction - parent_diff_value #(y=x-k)

  
summary_df = data.frame(input_file=input_file, sibling_mix=sibling_fraction,  
                        father_mix= father_fraction, mother_mix=mother_fraction, 
                        mendelian_error_rate = mendelian_error_rate)  
print(summary_df)
print(summary_path)
write_delim(summary_df, summary_path, delim='\t')






