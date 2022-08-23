suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(bbmle))
options(warn=-1)


option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input VAF sites at GroupA, B, C SNP loci", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"), 
  make_option(c("-r","--runmode"), type="character", default='all', help='[Optional] runmode, choose from single, joint, or all', metavar='character'),
  make_option(c("-u","--upd"), type="integer", default=1, help='[Optional] is set to 0 will not filter child variants with vaf=0 or 1 from a homo-ref + homo-alt parent loci. If set to 1, no filtering is applied', metavar='character'),
  make_option(c("-w","--wide"), type="integer", default=0, help='[Optional] is set to 1, produce final summary table in a wide format', metavar='character')
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


###############################################################################
##### vectorized likelihood calculations and negative-log-likelihood function 
# avoid log(0) = -Inf
check_float_min <-function(value){
  if(value==0){
    value = .Machine$double.xmin # to avoid falling in -Inf 
  }
  return(value)
}

###############################################################################
# simple mix of sibling only
likelihood_sibling<-function(alt, depth, ratio){
  p_00 = dbinom(alt, depth, prob=0)
  p_10 = dbinom(alt, depth, prob=ratio/2)
  p_01 = dbinom(alt, depth, prob= 0.5*(1-ratio))
  p_11 = dbinom(alt, depth, prob=0.5)
  prob_total = (p_00 + p_10 + p_01 + p_11 + 1e-100)/4
  prob_total = check_float_min(prob_total)
  return(log10(prob_total))
}

sibling_nll<-function(df, par){
  if(par>0){
    return(-sum(likelihood_sibling(df$alt, df$depth, par)))
  }else{
    return(10000000000000000)
  }
}

###############################################################################
# parent mix only
likelihood_parent_homohomo<-function(alt, depth, homo_parent, parent_diff){
  if(homo_parent=='F'){
    altfraction=min(1, max(0, (1-parent_diff)/2))
    p = dbinom(alt, depth, prob=altfraction)
    p = check_float_min(p)
    prob_total = (p + 1e-100)
    return(log10(prob_total))
  }else{
    altfraction=min(1, max(0, (1+parent_diff)/2))
    p = dbinom(alt, depth, prob=altfraction)
    p = check_float_min(p)
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

parent_diff_nll_single <-function(df, par){
  if(par[1] >= 0| par[1] <=1){
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
    p = check_float_min(p)

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
  p = check_float_min(p)

  prob_total = (p + 1e-100)
  return(log10(prob_total))
}

family_mix_nll<-function(df, k, x, z){
  if(x >= k & z >= 0 & (1-2*x+k-z >=0) & (z<=(1-2*x+k)/2)) {
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
# vectorized function for grid search
family_mix_nll_vect_par = Vectorize(family_mix_nll, vectorize.args = c('x', 'z'))


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
run_mode = opt$runmode
upd = opt$upd

# sanity check for input file existence and validity of runmode 
if(!(run_mode %in% c('all', 'single', 'joint'))){
  stop('--runmode should be either all, single, or joint. Exiting...')
}

if(!file.exists(input_file)){
  stop(paste0(input_file, ' is not found. Exiting...' ))
}

# read in count table
df <- read_delim(input_file, delim='\t', na=c('', 'NA'), col_types=cols(chrom=col_character(), pos=col_double(), alt=col_double(), depth=col_double(), vaf=col_double(), father_vaf=col_double(), mother_vaf=col_double()))

# filter sex chromosomes
df <- df %>% filter(!str_detect(chrom, pattern = 'chrX|X|chrY|Y')) # remove chrX for calculating MLE estimates

# filter out those regions where both parents are homoref/homoref for Mendelian error rate calculation
df_homorefhomoref <- df %>% filter(is.na(hetero_parent) & is.na(homoalt_parent))
denovo_error_rate = sum(df_homorefhomoref$alt)/ sum(df_homorefhomoref$depth) 

# identify homo alt homo ref
df_homoalthomoref <- df %>% filter(is.na(hetero_parent) & !is.na(homoalt_parent))


# filter out those regions where both parents are homoref/homoref
df_homoref_het <- df %>% filter(!is.na(hetero_parent) & is.na(homoalt_parent))


# Estimate MLE chimera ratio and confidence interval using mle function
# calcaulate using homo ref + homo_alt loci
df_homoalthomoref <- df %>% filter(is.na(hetero_parent) & !is.na(homoalt_parent))



if(upd==0){
  # if not interested in identifying upd, this filtering removes sites that are homozygous in the children (due to deletion)
  df_homoalthomoref = df_homoalthomoref %>% filter(vaf > 0 & vaf < 1) # added filter  to remove deletions
}

summary_df = list(input_file = input_file)

if(dim(df_homoalthomoref)[1] > 0 & dim(df_homoref_het)[1]>0){
  if(run_mode %in% c('all', 'joint')){
    parent_diff_result = mle2(parent_diff_nll, start=list(par=0), lower=list(par=-1), upper=list(par=1), data=list(df=df_homoalthomoref), method='Brent')
    parent_diff_value = coef(parent_diff_result)[['par']]

    # using the differnce between (mother-father=parent_diff) use this as a constraint, then find values for mother(x) and sibling(z) 
    mother_sibling_mle = mle2(family_mix_nll, start=list(x=max(parent_diff_value, 0), z=0), lower=list(x=max(parent_diff_value, 0), z=0), upper=list(x=(1+parent_diff_value)/2, z=0.5), fixed=list(k=parent_diff_value), data=list(df=df_homoref_het), method='L-BFGS-B')
    convergence_status = mother_sibling_mle@details$convergence
    # check for convergence of mle
    if(convergence_status!=0){
      print('Failed initial convergence of mle. Grid searching to estimate new initial values for optimizing')
      # use grid search to find a new starting point
      ll_space = as_tibble(expand.grid(x=0:100/100, z=0:50/100)) %>% filter(1-2*x+parent_diff_value-z >=0 & z<=(1-2*x+parent_diff_value)/2) %>% mutate(ll = family_mix_nll_vect_par(df=df_homoref_het, k=parent_diff_value, x=x, z=z))
      x_grid_estimate = as.double(ll_space[which(ll_space$ll == min(ll_space$ll)), 'x'])
      z_grid_estimate = as.double(ll_space[which(ll_space$ll == min(ll_space$ll)), 'z'])
      print(paste0('Initial guess for mother(x)=', round(x_grid_estimate, 4), ' sibling(z)=', round(z_grid_estimate, 4)))
      mother_sibling_mle = mle2(family_mix_nll, start=list(x=x_grid_estimate, z=z_grid_estimate), lower=list(x=max(parent_diff_value, 0), z=0), upper=list(x=(1+parent_diff_value)/2, z=0.5), fixed=list(k=parent_diff_value), data=list(df=df_homoref_het), method='L-BFGS-B', , control = list(maxit = 1e4, pgtol = 0, ndeps = c(1e-6, 1e-6), factr=0))
      convergence_status = mother_sibling_mle@details$convergence
    }


    mother_fraction = coef(mother_sibling_mle)[['x']]
    sibling_fraction = coef(mother_sibling_mle)[['z']]
    father_fraction = mother_fraction - parent_diff_value #(y=x-k)

    summary_df$child_contam_by_sibling_joint = sibling_fraction
    summary_df$child_contam_by_father_joint = father_fraction
    summary_df$child_contam_by_mother_joint = mother_fraction
    summary_df$convergence_joint=convergence_status
  }

  if(run_mode %in% c('all', 'single')){
     mother_single_mle = mle2(parent_diff_nll_single, start=list(par=0), lower=list(par=0), upper=list(par=1), data=list(df=df_homoalthomoref), method='Brent')
     father_single_mle = mle2(parent_diff_nll_single, start=list(par=0), lower=list(par=-1), upper=list(par=0), data=list(df=df_homoalthomoref), method='Brent')
     sibling_single_mle = mle2(sibling_nll, start=list(par=0), lower=list(par=0), upper=list(par=0.5), data=list(df=df_homoref_het), method='Brent')

     mother_fraction_single = coef(mother_single_mle)[['par']]
     father_fraction_single = -coef(father_single_mle)[['par']]
     sibling_fraction_single = coef(sibling_single_mle)[['par']]

     summary_df$child_contam_by_sibling = sibling_fraction_single
     summary_df$child_contam_by_father = father_fraction_single
     summary_df$child_contam_by_mother = mother_fraction_single

  }


    
  summary_df$denovo_error_rate = denovo_error_rate  



}
summary_df$groupA_father = df_homoalthomoref %>% filter(homoalt_parent=='F') %>% dim %>% `[[`(1)
summary_df$groupA_mother = df_homoalthomoref %>% filter(homoalt_parent=='M') %>% dim %>% `[[`(1)
summary_df$groupB_father = df_homoref_het %>% filter(hetero_parent=='F') %>% dim %>% `[[`(1)
summary_df$groupB_mother = df_homoref_het %>% filter(hetero_parent=='M') %>% dim %>% `[[`(1)

summary_df = data.frame(summary_df)

if(opt$wide==0){
  summary_df = summary_df %>% pivot_longer(names_to='type', values_to='value',-input_file) %>% select(-input_file)
}

write_delim(summary_df, summary_path, delim='\t')






