# infer the most likely chimeric ratio from data
# 2019.08.07 cjyoon
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))


option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input VAF sites at HET/HOMREF parental sites", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"),
  make_option(c("-s","--seq"), type="character", default="0,0.5,0.001", help="[Optional] seq function in R for chimera ratio search space", metavar="character"),
  make_option(c("-f","--father_hetero_region"), type="character", default=NULL, help="Father hetero Region to fit MLE", metavar="character"),
  make_option(c("-m","--mother_hetero_region"), type="character", default=NULL, help="mother hetero Region to fit MLE", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

likelihood<-function(alt, depth, ratio){
  p_10 = dbinom(as.integer(alt), as.integer(depth), prob=ratio/2)
  p_01 = dbinom(as.integer(alt), as.integer(depth), prob=0.5*(1-ratio))
  prob_total = (p_10 + p_01 + 1e-100)/2
  return(log10(prob_total))
}

chimeria_mle_estimate<-function(alt_count, depth, seq_start, seq_end, seq_by, title=NA, output_path=NA){
  
  likelihood_space=rep(NA, 501)
  estimated_ratio_space = seq(seq_start, seq_end, by=seq_by)
  i=1
  result <- sapply(estimated_ratio_space, function(estimated_ratio){
    args_list = list(alt_count = alt_count, depth = depth, ratio=estimated_ratio)
    
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
  print(mle_estimate)
  p<-ggplot(df, aes(x=estimated_ratio_space, y=loglikelihood)) +
    geom_point() + geom_line() + 
    geom_vline(xintercept = mle_estimate, linetype='dashed', color='red', size=1) + 
    annotate('text', x=mle_estimate + 0.02, y=min(df$loglikelihood), label=paste0('MLE estimate = ', as.character(mle_estimate))) + 
    xlab('Estimated chimeric ratio') + ylab('Log Likelihood') + 
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

seq_string = opt$seq

seq_string= as.double(str_split(string=seq_string, pattern=',')[[1]])
seq_start= seq_string[1]
seq_end = seq_string[2]
seq_by = seq_string[3]

father_region =  opt$father_hetero_region #'~/Projects/chimera/chimera_triplet/father_discordant_regions.bed'
father_region_df = read_delim(father_region, delim='\t', col_names=c('chrom', 'start', 'end'))
mother_region = opt$mother_hetero_region #'~/Projects/chimera/chimera_triplet/mother_discordant_regions.bed'#
mother_region_df = read_delim(mother_region, delim='\t', col_names=c('chrom', 'start', 'end'))



output_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file),'_discordant_regions.mle.pdf')))
summary_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file), '_discordant_regions.summary.tsv')))


df <- read_delim(input_file, delim='\t')

print('filter X')
df <- df %>% filter(!str_detect(chrom, pattern = 'chrX|X')) # remove chrX for calculating MLE estimates

father_df = df %>% filter(hetero_parent=='F')
mother_df = df %>% filter(hetero_parent == 'M')

father_df_roi = data.frame()
for(region_i in 1:dim(father_region_df)[1]){
  chromosome_i=father_region_df$chrom[region_i]
  start_i = as.double(father_region_df$start[region_i])
  end_i = as.double(father_region_df$end[region_i])
  father_df_roi_i = father_df %>% filter(chrom==chromosome_i) %>% filter(pos>start_i) %>% filter(pos < end_i)
  father_df_roi = rbind(father_df_roi, father_df_roi_i)
}

mother_df_roi = data.frame()
for(region_i in 1:dim(mother_region_df)[1]){
  chromosome_i=mother_region_df$chrom[region_i]
  start_i = as.double(mother_region_df$start[region_i])
  end_i = as.double(mother_region_df$end[region_i])
  mother_df_roi_i = mother_df %>% filter(chrom==chromosome_i) %>% filter(between(pos, start_i, end_i))
  mother_df_roi = rbind(mother_df_roi, mother_df_roi_i)
}

df = rbind(father_df_roi, mother_df_roi)
print(df)

# Estimate MLE from optim
print('Using optim function to estimate chimeric ratio')
loss_function = function(df, par){
  if(par[1] > 0){
    # probability is always a positive number
    return(-sum(likelihood(df$alt, df$depth, par)))
  }else{
    # tell optim to stay away as probability can only be a positive number
    return(10000000000000000)
  }
}

result <- optim(par=0.25, d=df, fn=loss_function, method='Brent', lower=0, upper=0.5, hessian = T)
hessian <- result$hessian
hessian.inv <- solve(hessian)
parameter.se <- sqrt(diag(hessian.inv))


mle <- result$par
upper_ci <- result$par - 1.96 * parameter.se
lower_ci <- result$par + 1.96 * parameter.se

summary_df = data.frame(input_file=input_file, mle=mle, se = parameter.se, upper_ci = upper_ci, lower_ci = lower_ci)  
write_delim(summary_df, summary_path, delim='\t')
print(summary_df)


# # Estimate MLE chimera ratio
# system.time(result <-chimeria_mle_estimate(df$alt, df$depth, seq_start, seq_end, seq_by, title=basename(output_path), output_path=output_path))

# # Plot histogram of VAFs
# output_histogram = normalizePath(file.path(opt$output_dir, paste0(basename(output_path), '.hist.pdf')))
# ggplot(df, aes(x=vaf)) + 
#   geom_histogram(aes(y=..density..), colour="black", fill="white", breaks=seq(0, 1, by = 0.01))+ 
#   geom_density(fill='skyblue', alpha=0.3) + 
#   xlab('VAF') + ylab('Density') + 
#   ggtitle(basename(input_file))
# ggsave(output_histogram)




