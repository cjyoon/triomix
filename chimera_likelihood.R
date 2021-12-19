# infer the most likely chimeric ratio from data
# 2019.08.07 cjyoon
# 2021.12.19 cjyoon added/fixed parent mix
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
library(optimParallel)
options(warn=-1)


option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input VAF sites at HET/HOMREF parental sites", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"),
  make_option(c("-r","--run_mode"), type="character", default='plot', help="Run mode, if set to 'plot', then full calculation of log likelihood in the parameter space, if set to 'optim', then calculate quickly using optim with confidence interval calculation, if set to 'all', then run both full calculation and quick calculation with optim", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


####### LIKELIHOOD FUNCTIONS
likelihood<-function(alt, depth, ratio){
  p_00 = dbinom(alt, depth, prob=0)
  p_10 = dbinom(alt, depth, prob=ratio/2)
  p_01 = dbinom(alt, depth, prob=0.5*(1-ratio))
  p_11 = dbinom(alt, depth, prob=0.5)
  prob_total = (p_00 + p_10 + p_01 + p_11 + 1e-100)/4
  return(log10(prob_total))
}
likelihood_parentmix<-function(alt, depth, ratio){
  # p_10 = dbinom(alt, depth, prob=ratio/2)
  p_01 = dbinom(alt, depth, prob=0.5*(1-ratio))
  prob_total = (p_01 + 1e-100)/2
  return(log10(prob_total))
}

likelihood_homoalt_homoref<-function(alt, depth, ratio){
  p = dbinom(alt, depth, prob=0.5*(1+ratio))
  prob_total = (p + 1e-100)/4
  return(log10(prob_total))
}
likelihood_homoalt_homoref_combined<-function(alt, depth, homoalt_parent, parent_diff){
  # parent_diff = (mother % - father %)
  if(homoalt_parent=='F'){
    altfraction=(1-parent_diff)/2
  }else if(homoalt_parent=='M'){
    altfraction=(1+parent_diff)/2
  }
  p = dbinom(alt, depth, prob=altfraction)
  prob_total = (p + 1e-100)
  return(log10(prob_total))
}
likelihood_homoalt_homoref_combined_vec = Vectorize(likelihood_homoalt_homoref_combined)


likelihood_het_homoref_parent_combined<-function(alt, depth, hetalt_parent, parent_diff, mother_fraction){
  if(hetalt_parent=='F'){
    p1 = dbinom(alt, depth, prob=(1-mother_fraction)/2)
    p0 = dbinom(alt, depth, prob=(mother_fraction - parent_diff)/2)
  }else if(hetalt_parent=='M'){
    p1 = dbinom(alt, depth, prob=(1-mother_fraction + parent_diff)/2)
    p0 = dbinom(alt, depth, prob=(mother_fraction)/2)
  }
  prob_total = (p0 + p1 + 1e-100)/2
  return(log10(prob_total))
}
likelihood_het_homoref_parent_combined_vec = Vectorize(likelihood_het_homoref_parent_combined)

###### LOSS FUNCTIONS FOR optim
loss_function = function(df, par){
  if(par[1] > 0){
    # probability is always a positive number
    return(-sum(likelihood(df$alt, df$depth, par)))
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}

loss_function_parent = function(df, par){
  if(par[1] > 0){
    # probability is always a positive number
    return(-sum(likelihood_parentmix(df$alt, df$depth, par)))
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}

loss_function_parent_homoalt_homoref = function(df, par){
  if(par[1] > 0){
    # probability is always a positive number
    return(-sum(likelihood_homoalt_homoref(df$alt, df$depth, par)))
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}

loss_function_parent_homoalt_homoref_combined = function(df, par){
  if(par[1] > 0){
    # probability is always a positive number
    return(-sum(likelihood_homoalt_homoref_combined_vec(df$alt, df$depth, df$homoalt_parent, par)))
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}

loss_function_parent_het_homoref_combined = function(df, par){
  if(par[1] > 0){
    # probability is always a positive number
    return(-sum(likelihood_het_homoref_parent_combined_vec(df$alt, df$depth, df$hetero_parent, df$parent_diff, par)))
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}

chimeria_mle_estimate<-function(alt_count, depth, title=NA, output_path=NA){
  
  estimated_ratio_space = seq(0, 0.5, by=0.001)
  i=1
  result <- sapply(estimated_ratio_space, function(estimated_ratio){
    args_list = list(alt = alt_count, depth = depth, ratio=estimated_ratio)
    
    log_likelihood_list = unlist(pmap(args_list, likelihood))
    
    sum_likelihood = sum(log_likelihood_list)
    sum_likelihood
  })
  df <- data.frame(chimera_ratio = estimated_ratio_space, loglikelihood=result)
  mle_estimate = (estimated_ratio_space[which(max(result)==result)])
  # if(!is.na(output_path)){
  #   pdf(output_path)
  # }
  if(is.na(title)) title=''
  
  p<-ggplot(df, aes(x=estimated_ratio_space, y=loglikelihood)) +
    geom_point() + geom_line() + 
    geom_vline(xintercept = mle_estimate, linetype='dashed', color='red', size=1) + 
    annotate('text', x=mle_estimate + 0.07, y=min(df$loglikelihood), label=paste0('MLE estimate = ', as.character(mle_estimate))) + 
    xlab('Estimated chimeric ratio') + ylab('Log Likelihood') +
    theme_classic(base_size = 24) + 
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    ggtitle(title)
  ggsave(output_path)
  # if(!is.na(output_path)){ 
  #   dev.off()
  # }
  result_path = gsub('.pdf', '.txt', output_path)
  write_delim(df, paste0(output_path,  '.data.tsv'), delim='\t')
  writeLines(paste0(input_file, '\t', as.character(mle_estimate)), result_path)
  
  p
}

######################################
# Read in input file and specify output files

if(opt$output_dir=='.' || is.na(opt$output_dir)){
  opt$output_dir = getwd() 
}

input_file = opt$input_file
output_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.mle.pdf')))
summary_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.summary.tsv')))
run_mode = opt$run_mode
df <- read_delim(input_file, delim='\t', na=c('', 'NA'), col_types=cols(pos=col_double(), alt=col_double(), depth=col_double(), vaf=col_double(), father_vaf=col_double(), mother_vaf=col_double()))
df <- df %>% filter(!str_detect(chrom, pattern = 'chrX|X')) # remove chrX for calculating MLE estimates

# filter out those regions where both parents are homoref/homoref
df_homorefhomoref <- df %>% filter(is.na(hetero_parent) & is.na(homoalt_parent))
#calculate error rate
mendelian_error_rate = sum(df_homorefhomoref$alt)/ sum(df_homorefhomoref$depth) * 2 # multiplied by 2 since diploid genome
# identify homo alt homo ref
df_homoalthomoref <- df %>% filter(is.na(hetero_parent) & !is.na(homoalt_parent))


# filter out those regions where both parents are homoref/homoref
df_homoref_het <- df %>% filter(!is.na(hetero_parent) & is.na(homoalt_parent))


# Estimate MLE chimera ratio and confidence interval using optim function
if(run_mode %in% c('optim', 'all')){
  
  print('Using optim function to estimate chimeric ratio')
  result <- optim(par=0.25, d=df_homoref_het, fn=loss_function, method='Brent', lower=0, upper=0.5, hessian = T)
  hessian <- result$hessian
  hessian.inv <- solve(hessian)
  parameter.se <- sqrt(diag(hessian.inv))
  
  
  # mother_variants_df = df_homoref_het %>% filter(hetero_parent=='M')
  # father_variants_df = df_homoref_het %>% filter(hetero_parent=='F')
  
  # mother_het_result = optim(par=0.5, d=mother_variants_df, fn=loss_function_parent, method='Brent', lower=0, upper=1, hessian = T)
  # father_het_result = optim(par=0.5, d=father_variants_df, fn=loss_function_parent, method='Brent', lower=0, upper=1, hessian = T)
  
  
  # calcaulate using homo ref + homo_alt loci
  df_homoalthomoref <- df %>% filter(is.na(hetero_parent) & !is.na(homoalt_parent))
  mother_homoalt = df_homoalthomoref %>% filter(homoalt_parent=='M')
  father_homoalt = df_homoalthomoref %>% filter(homoalt_parent=='F')
  
  # mother_homo_result = optim(par=0, d=mother_homoalt, fn=loss_function_parent_homoalt_homoref_combined, method='Brent', lower=-1, upper=1, hessian = T)
  # father_homo_result = optim(par=0, d=father_homoalt, fn=loss_function_parent_homoalt_homoref_combined, method='Brent', lower=-1, upper=1, hessian = T)
  parent_diff = optim(par=0, d=df_homoalthomoref, fn=loss_function_parent_homoalt_homoref_combined, method='Brent', lower=-1, upper=1, hessian = T)

  
  df_homoref_het$parent_diff = parent_diff$par
  mother_fraction_result = optim(par=0, d=df_homoref_het, fn=loss_function_parent_het_homoref_combined, method='Brent', lower=0, upper=1, hessian=T)
  mother_fraction = mother_fraction_result$par
  father_fraction = mother_fraction - parent_diff$par
  
  mle <- result$par
  upper_ci <- result$par - 1.96 * parameter.se
  lower_ci <- result$par + 1.96 * parameter.se
  
  summary_df = data.frame(input_file=input_file, mle=mle, se = parameter.se, upper_ci = upper_ci, lower_ci = lower_ci,
                          # mother_het_result=mother_het_result$par, father_het_result=father_het_result$par, 
                          mother_mix= mother_fraction, father_mix=father_fraction, 
                          mendelian_error_rate = mendelian_error_rate)  
  print(summary_df)
  print(summary_path)
  write_delim(summary_df, summary_path, delim='\t')
}


# Estimate MLE chimera ratio in full calculation mode
if(run_mode %in% c('plot', 'all')){
  print('Plot mode, calculating log likelihood between 0~0.5 chimeric ratio space')
  system.time(result <-chimeria_mle_estimate(df_homoref_het$alt, df_homoref_het$depth, title=input_file, output_path=output_path))
  
  # Plot histogram of VAFs
  output_histogram = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.hist.pdf')))
  ggplot(df, aes(x=vaf)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", breaks=seq(0, 1, by = 0.01))+ 
    geom_density(fill='skyblue', alpha=0.3) + 
    xlab('VAF') + ylab('Density') + 
    ggtitle(basename(input_file))
  ggsave(output_histogram)
}






