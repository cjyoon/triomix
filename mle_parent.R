suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(bbmle))
options(warn=-1)


option_list = list(
  make_option(c("-i","--input_parent_counts"), type="character", default=NULL, help="Readcounts at autosomal SNPs where child is homo-alt", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"), 
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


likelihood <- function(alt, depth, ratio){
  p = check_float_min(dbinom(alt, depth, prob=ratio))
  return(log10(p))
}
likelihood_vec = Vectorize(likelihood, vectorize.args = c('alt', 'depth'))

father_contam_by_child_nll <-function(df, par){
  if(par>=0 & par <=1){
    nll_sum = -sum(likelihood_vec(df$father_alt, df$father_depth, par))
    return(nll_sum)
  }else{
    return(10000000000000000)
  }
}

mother_contam_by_child_nll <-function(df, par){
  if(par>=0 & par <=1){
    nll_sum = -sum(likelihood_vec(df$mother_alt, df$mother_depth, par))
    return(nll_sum)
  }else{
    return(10000000000000000)
  }
}


contam_in_father_nll <- function(df, par){
  if(par>=0 & par <=1){
    nll_sum = -sum(likelihood_vec(df$father_alt, df$father_depth, par))
    return(nll_sum)
  }else{
    return(10000000000000000)
  }
}


contam_in_mother_nll <-function(df, par){
  if(par>=0 & par <=1){
    nll_sum = -sum(likelihood_vec(df$mother_alt, df$mother_depth, par))
    return(nll_sum)
  }else{
    return(10000000000000000)
  }
}


df = read_delim(opt$input_parent_counts, delim='\t', 
                na=c('', 'NA'), 
                col_types=cols(chrom=col_character(), 
                               pos=col_double(), 
                               alt=col_double(), 
                               depth=col_double(), 
                               vaf=col_double(), 
                               father_vaf=col_double(), 
                               mother_vaf=col_double()))

# filter out sex chromosome for analysis
df = df %>% filter(!str_detect(chrom, '^X$|^chrX$|^Y$|^chrY$'))

##############################################
# Filter df for GroupD variants
df_mother_by_child = df %>% filter(mother_vaf < 1) %>% filter(mother_depth - mother_alt >1)
df_father_by_child = df %>% filter(father_vaf < 1) %>% filter(father_depth - father_alt >1)


# MLE fitting 
if(dim(df_mother_by_child)[1]>0){
  mother_by_child = mle2(contam_in_mother_nll, start=list(par=0.5), lower=list(par=0.5), upper=list(par=1), data=list(df=df_mother_by_child), method='Brent')
  mother_by_child_res = coef(mother_by_child)[['par']]
}else{
  mother_by_child_res = NA # if not enough data, return NA. This occurs if there is too much contmination in mother
}


if(dim(df_father_by_child)[1]>0){
  father_by_child = mle2(contam_in_father_nll, start=list(par=0.5), lower=list(par=0.5), upper=list(par=1), data=list(df=df_father_by_child), method='Brent')
  father_by_child_res = coef(father_by_child)[['par']]
}else{
  father_by_child_res = NA # if not enough data, return NA. This occurs if there is too much contmination in mother
}


w_mother  = 2*mother_by_child_res -1
w_father  = 2*father_by_child_res -1

###################################################################
# Filter df for GroupE variants
df_father_by_mother =  df %>% filter(father_vaf < 1) %>% filter(father_depth - father_alt >1) %>% filter(mother_vaf==1)
df_mother_by_father = df %>% filter(mother_vaf < 1) %>% filter(mother_depth - mother_alt >1) %>% filter(father_vaf==1)




if(dim(df_mother_by_father)[1]>0){
  mother_by_father = mle2(contam_in_mother_nll, start=list(par=0.5), lower=list(par=0.5), upper=list(par=1), data=list(df=df_mother_by_father), method='Brent')
  mother_by_father_res = coef(mother_by_father)[['par']]
}else{
  mother_by_father_res=NA
}

if(dim(df_father_by_mother)[1]>0){
  father_by_mother = mle2(contam_in_father_nll, start=list(par=0.5), lower=list(par=0.5), upper=list(par=1), data=list(df=df_father_by_mother), method='Brent')
  father_by_mother_res = coef(father_by_mother)[['par']]
}else{
  father_by_mother_res=NA
}



v_mother  = 2*mother_by_father_res -1
v_father  = 2*father_by_mother_res -1


# write output file
output_path = normalizePath(file.path(opt$output_dir, paste0(basename(opt$input_parent_counts), '.summary.tsv')))

output_df = data.frame( input_file = opt$input_parent_counts, 
                        mother_contam_by_child = w_mother,
                        father_contam_by_child = w_father,
                        mother_contam_by_father = v_mother, 
                        father_contam_by_mother = v_father, 
                        groupD_mother = dim(df_mother_by_child)[1], 
                        groupD_father = dim(df_father_by_child)[1], 
                        groupE_mother = dim(df_mother_by_father)[1], 
                        groupE_father = dim(df_father_by_mother)[1])

if(opt$wide==0){
  output_df = output_df %>% pivot_longer(names_to='type', values_to='value',-input_file) %>% select(-input_file)
}

write_delim(output_df, output_path, delim='\t')



