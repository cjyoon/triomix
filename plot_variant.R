# plot VAF, depth etc for triomix
# 2022.01.01 cjyoon draw whole genome wide vaf plot

suppressMessages(library(tidyverse))
suppressMessages(library(grid))
suppressMessages(library(optparse))
options(warn=-1)


option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input counts file", metavar="character"),
  make_option(c("-r","--reference"), type="character", default=NULL, help="Reference FASTA file. Must be indexed", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"),
  make_option(c("-s","--subsample"), type="double", default=1, help="Subsampling ratio", metavar="character")

)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


counts_path = opt$input_file
reference_path = opt$reference
fai_path = paste0(reference_path, '.fai')
subsample_ratio = opt$subsample
print(subsample_ratio)
output_path = normalizePath(file.path(opt$output_dir, paste0(basename(counts_path), '.plot.pdf')))



counts_df = read_delim(counts_path, delim='\t')
reference_fai = read_delim(fai_path, delim='\t', col_names = c('name', 'length', 'offset', 'linebases', 'linewidth', 'qualoffset'))
reference_fai = reference_fai %>% filter(str_detect(name, '^chr[0-9XY]+$|^[0-9XY]+$'))
reference_fai = reference_fai %>% mutate(cumpos = cumsum(length), lag1=lag(cumpos), lag1=replace_na(lag1, 0))
reference_fai_dict = setNames(reference_fai$lag1, reference_fai$name)
total_genome_length = sum(reference_fai$length)

chrom_to_numeric <- function(chromosome, position, reference_fai_dict){
  # change chromosomal position into a single numeric value for plotting
  return(reference_fai_dict[[chromosome]] + position)
}
chrom_to_numeric_vec = Vectorize(chrom_to_numeric, vectorize.args = c('chromosome', 'position'))


counts_df_numpos = counts_df %>% mutate(numeric_pos = chrom_to_numeric_vec(chrom, pos, reference_fai_dict))



plot_vaf <- function(count_df, parent, reference_fai_dict, total_genome_length, subsample_ratio=1, plot_label=''){
  parent_color = c('father'='red', 'mother'='blue')
  color = parent_color[[parent]]
  pushViewport(viewport(x=0.05, y=0.2, width=0.9, height=0.7, just=c('left', 'bottom')))
  grid.rect()
  grid.text(plot_label, x = 0, y=1, just=c('left', 'bottom'))
  for(i in 2:length(reference_fai_dict)){
    fai_length = reference_fai_dict[i]
    prev_length = reference_fai_dict[i-1]
    chrom_text_pos = (fai_length + prev_length)/2
    
    grid.lines(x = c(fai_length/total_genome_length, fai_length/total_genome_length), y=c(0, 1))
    grid.text(label =names(reference_fai_dict)[i-1], x= chrom_text_pos/total_genome_length, y=-0.05, rot = 45, just=c('right', 'center'))
    
  }
  count_df = count_df %>% sample_frac(subsample_ratio)
  
  grid.yaxis(at=c(0, 0.5, 1), label=c(0, 0.5, 1))
  grid.text('VAF', rot = 90, x=-0.04, y=0.5, just=c('center', 'bottom'))
  grid.points(
    x=count_df$numeric_pos/total_genome_length, 
    y=count_df$vaf, 
    gp=gpar(col=color, alpha=0.3), pch=16, size = unit(0.5, 'char')
  )
  popViewport(1)
}

plot_depth <- function(count_df, reference_fai_dict, total_genome_length, subsample_ratio=1){
  pushViewport(viewport(x=0.05, y=0.2, width=0.9, height=0.7, just=c('left', 'bottom')))
  grid.rect()
  grid.text('Depth', x = 0, y=1, just=c('left', 'bottom'))
  for(i in 2:length(reference_fai_dict)){
    fai_length = reference_fai_dict[i]
    prev_length = reference_fai_dict[i-1]
    chrom_text_pos = (fai_length + prev_length)/2
    
    grid.lines(x = c(fai_length/total_genome_length, fai_length/total_genome_length), y=c(0, 1))
    grid.text(label =names(reference_fai_dict)[i-1], x= chrom_text_pos/total_genome_length, y=-0.05, rot = 45, just=c('right', 'center'))
  }
  count_df = count_df %>% sample_frac(subsample_ratio)
  mean_depth = mean(count_df$depth, na.rm=T)
  
  grid.points(
    x=count_df$numeric_pos/total_genome_length, 
    y=count_df$depth/(mean_depth*2), 
    gp=gpar(alpha=0.3), pch=16,  size = unit(0.5, 'char')
  )
  grid.yaxis(at=seq(0, mean_depth*2, 10)/(mean_depth*2), label=seq(0, mean_depth*2, 10))
  grid.text('Depth', rot = 90, x=-0.04, y=0.5, just=c('center', 'bottom'))
  
  popViewport(1)
}




# select homo-ref(mother)/homo-alt(father)
homoref_homo_alt_father = counts_df_numpos %>% filter(father_vaf==1 & mother_vaf==0)

# select homo-ref(father)/homo-alt(mother)
homoref_homo_alt_mother = counts_df_numpos %>% filter(mother_vaf==1 & father_vaf==0)


# select homo-ref(mother)/het(father)
homoref_het_father = counts_df_numpos %>% filter(mother_vaf==0 & father_vaf > 0.4 & father_vaf < 0.6)

# select homo-ref(father)/het(mother)
homoref_het_mother = counts_df_numpos %>% filter(father_vaf==0 & mother_vaf > 0.4 & mother_vaf < 0.6)


grid.newpage()
pdf(output_path, width=20, height=10)
pushViewport(viewport(x=0,  y=0, width=1, height=0.2, just=c('left', 'bottom')))
plot_vaf(homoref_het_father, 'father', reference_fai_dict, total_genome_length=total_genome_length, subsample_ratio=subsample_ratio, plot_label='homo-ref (mother) + het (father)')
popViewport(1)
pushViewport(viewport(x=0,  y=0.20, width=1, height=0.2, just=c('left', 'bottom')))
plot_vaf(homoref_het_mother, 'mother', reference_fai_dict, total_genome_length=total_genome_length, subsample_ratio=subsample_ratio, plot_label='homo-ref (father) + het (mother)')
popViewport(1)

pushViewport(viewport(x=0,  y=0.40, width=1, height=0.20, just=c('left', 'bottom')))
plot_vaf(homoref_homo_alt_father, 'father', reference_fai_dict, total_genome_length=total_genome_length, subsample_ratio=subsample_ratio, plot_label='homo-ref (mother) + homo-alt (father)')
popViewport(1)

pushViewport(viewport(x=0,  y=0.60, width=1, height=0.2, just=c('left', 'bottom')))
plot_vaf(homoref_homo_alt_mother, 'mother', reference_fai_dict, total_genome_length=total_genome_length, subsample_ratio=subsample_ratio, plot_label='homo-ref (father) + homo-alt (mother)')
popViewport(1)


pushViewport(viewport(x=0,  y=0.80, width=1, height=0.2, just=c('left', 'bottom')))
plot_depth(bind_rows(homoref_het_father, homoref_het_mother, homoref_homo_alt_father, homoref_homo_alt_mother), reference_fai_dict, total_genome_length=total_genome_length, subsample_ratio=subsample_ratio)
popViewport(1)
dev.off()


