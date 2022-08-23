
suppressMessages(library(PSCBS))
suppressMessages(library(DNAcopy))
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
options(warn=-1)





option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input counts file", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"),
  make_option(c("-s","--segment_length"), type="double", default=2, help="minimum segment length in mega bases(MB)", metavar="integer")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


counts_path = opt$input_file
counts_df = read_delim(counts_path, delim='\t')
segment_length = opt$segment_length



convert_chr_to_integer <- function(x){
  # written since PSCBS requires chr to be in integer
  if(str_detect(x, '[0-9]+')){
    as.integer(gsub(pattern = 'chr', replacement = '', x=x))
  }else{
    23 # chrX -> 23
  }
}

convert_chr_to_integer_vec = Vectorize(convert_chr_to_integer)

convert_integer_to_chrom <- function(x){
  if(x==23){
    return('chrX')
  }else{
    return(paste0('chr', x))
  }
}
convert_integer_to_chrom_vec = Vectorize(convert_integer_to_chrom)







counts_df = counts_df %>% mutate(chrom=convert_chr_to_integer_vec(chrom))
counts_df = counts_df %>% rename(chromosome=chrom, x=pos, y=vaf)

df_homoalt_father = counts_df %>% filter(father_vaf==1 & mother_vaf==0) %>% select(chromosome, x, y)
df_homoalt_mother = counts_df %>% filter(mother_vaf==1 & father_vaf==0) %>% select(chromosome, x, y)

print(dim(df_homoalt_father))
print(dim(df_homoalt_mother))

homorefalt_segments = data.frame()
for(chrom in unique(counts_df$chromosome)){
  
  df_homoalt = list(father=df_homoalt_father %>% filter(chromosome==chrom),
                    mother=df_homoalt_mother %>% filter(chromosome==chrom))
  
  for(parent in c('father', 'mother')){
    df_homoalt_parent = df_homoalt[[parent]]
    if(dim(df_homoalt_parent)[1]>1){ # run only if more than 1 data point exists
      df_homoalt_parent = dropSegmentationOutliers(df_homoalt_parent)
    
    gaps <- findLargeGaps(df_homoalt_parent, minLength=segment_length * 1e+06)
    knownSegments <- gapsToSegments(gaps)
    
    fit <- segmentByCBS(df_homoalt_parent, knownSegments = knownSegments, seed = 48879, verbose = -10)
    segmented = getSegments(fit, simplify = TRUE)
    
    segmented = segmented %>% select(-sampleName)
    segmented$parent = parent
    homorefalt_segments = bind_rows(homorefalt_segments, segmented)
    }
  }
}

if(dim(homorefalt_segments)[1]>0){
  homorefalt_segments = homorefalt_segments %>% mutate(chromosome=convert_integer_to_chrom_vec(chromosome)) %>% arrange(desc(parent))
}else{
  print('segmentation not performed due to lack of GroupA variants')
}


output_path = normalizePath(file.path(opt$output_dir, paste0(basename(counts_path), '.upd.segments.tsv')))



write_delim(homorefalt_segments, output_path, delim='\t')


