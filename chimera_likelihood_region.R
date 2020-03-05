# infer the most likely chimeric ratio from data
# 2019.08.07 cjyoon
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))



option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input VAF sites at HET/HOMREF parental sites", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"),
  make_option(c("-s","--seq"), type="character", default="0,0.5,0.001", help="[Optional] seq function in R for chimera ratio search space", metavar="character"),
  make_option(c("-r","--region"), type="character", default=NULL, help="Region to fit MLE", metavar="character"),
  make_option(c("-p","--hetero_parent"), type="character", default=NULL, help="Hetero parent (F=father, M=mother)", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

likelihood<-function(alt_count, depth, ratio){
  p_10 = dbinom(alt_count, depth, prob=ratio/2)
  p_01 = dbinom(alt_count, depth, prob=0.5*(1-ratio))
  prob_total = (p_10 + p_01 + 1e-100)/4
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

region = opt$region
chromosome = str_split(region, ':')[[1]][1]
start_end = str_split(region, ':')[[1]][2]
region_start = as.integer(str_split(start_end, '-')[[1]][1])
region_end = as.integer(str_split(start_end, '-')[[1]][2])

hp = opt$hetero_parent

region_string = gsub(':', '_', region)
output_path = normalizePath(file.path(opt$output_dir, paste0(basename(input_file),'_', region_string, '_', hp,  '.mle.pdf')))


df <- read_delim(input_file, delim='\t')

print('filter X')
df <- df %>% filter(!str_detect(chrom, pattern = 'chrX|X')) # remove chrX for calculating MLE estimates
print(df)
# filter by region
print('region filter')
df <- df %>% filter(chrom==chromosome) %>% filter(between(pos, region_start, region_end))
print(df)
# filter by hetero parent
print('filter heteroparent')
df <- df %>% filter(hetero_parent == hp)

print(dim(df))
# Estimate MLE chimera ratio
system.time(result <-chimeria_mle_estimate(df$alt, df$depth, seq_start, seq_end, seq_by, title=basename(output_path), output_path=output_path))

# Plot histogram of VAFs
output_histogram = normalizePath(file.path(opt$output_dir, paste0(basename(output_path), '.hist.pdf')))
ggplot(df, aes(x=vaf)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", breaks=seq(0, 1, by = 0.01))+ 
  geom_density(fill='skyblue', alpha=0.3) + 
  xlab('VAF') + ylab('Density') + 
  ggtitle(basename(input_file))
ggsave(output_histogram)




