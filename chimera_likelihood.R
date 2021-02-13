# infer the most likely chimeric ratio from data
# 2019.08.07 cjyoon
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))



option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input VAF sites at HET/HOMREF parental sites", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"),
  make_option(c("-r","--run_mode"), type="character", default='plot', help="Run mode, if set to 'plot', then full calculation of log likelihood in the parameter space, if set to 'optim', then calculate quickly using optim with confidence interval calculation, if set to 'all', then run both full calculation and quick calculation with optim )", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

likelihood<-function(alt, depth, ratio){
  p_00 = dbinom(alt, depth, prob=0)
  p_10 = dbinom(alt, depth, prob=ratio/2)
  p_01 = dbinom(alt, depth, prob=0.5*(1-ratio))
  p_11 = dbinom(alt, depth, prob=0.5)
  prob_total = (p_00 + p_10 + p_01 + p_11 + 1e-100)/4
  return(log10(prob_total))
}

chimeria_mle_estimate<-function(alt_count, depth, title=NA, output_path=NA){
  
  likelihood_space=rep(NA, 501)
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
input_file = opt$input_file
output_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.mle.pdf')))
summary_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.summary.tsv')))
run_mode = opt$run_mode
df <- read_delim(input_file, delim='\t')
df <- df %>% filter(!str_detect(chrom, pattern = 'chrX|X')) # remove chrX for calculating MLE estimates

# Estimate MLE chimera ratio and confidence interval using optim function
if(run_mode %in% c('optim', 'all')){
  
  loss_function = function(df, par){
    if(par[1] > 0){
      # probability is always a positive number
      return(-sum(likelihood(df$alt, df$depth, par)))
    }else{
      # tell optim to stay away as probability can only be a positive number
      return(10000000000000000)
    }
  }
  print('Using optim function to estimate chimeric ratio')
  result <- optim(par=0.25, d=df, fn=loss_function, method='Brent', lower=0, upper=0.5, hessian = T)
  hessian <- result$hessian
  hessian.inv <- solve(hessian)
  parameter.se <- sqrt(diag(hessian.inv))
  
  
  mle <- result$par
  upper_ci <- result$par - 1.96 * parameter.se
  lower_ci <- result$par + 1.96 * parameter.se
  
  summary_df = data.frame(input_file=input_file, mle=mle, se = parameter.se, upper_ci = upper_ci, lower_ci = lower_ci)  
  write_delim(summary_df, summary_path, delim='\t')
}


# Estimate MLE chimera ratio in full calculation mode
if(run_mode %in% c('plot', 'all')){
  print('Plot mode, calculating log likelihood between 0~0.5 chimeric ratio space')
  system.time(result <-chimeria_mle_estimate(df$alt, df$depth, title=input_file, output_path=output_path))
  
  # Plot histogram of VAFs
  output_histogram = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '.hist.pdf')))
  ggplot(df, aes(x=vaf)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", breaks=seq(0, 1, by = 0.01))+ 
    geom_density(fill='skyblue', alpha=0.3) + 
    xlab('VAF') + ylab('Density') + 
    ggtitle(basename(input_file))
  ggsave(output_histogram)
}






