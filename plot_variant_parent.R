# plot VAF, depth etc for triomix for the contamination in the parents
suppressMessages(library(tidyverse))
suppressMessages(library(grid))
suppressMessages(library(optparse))
options(warn=-1)


option_list = list(
  make_option(c("-i","--input_file"), type="character", default=NULL, help="Input counts file", metavar="character"),
  make_option(c("-r","--reference"), type="character", default=NULL, help="Reference FASTA file. Must be indexed", metavar="character"),
  make_option(c("-o","--output_dir"), type="character", default="./", help="[Optional] Output directory (default: ./)", metavar="character"),
  make_option(c("-d","--downsample"), type="double", default=1, help="Downsampling ratio", metavar="character"),
  make_option(c("-a","--autosome"), type="logical", default=F, action='store_true', help="Draw autosome only", metavar="character"),
  make_option(c("-f","--output_format"), type="character", default='pdf', help="Output format. Options: pdf, png, jpg. Default=pdf", metavar="character")

  
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)



counts_path = opt$input_file
reference_path = opt$reference
fai_path = paste0(reference_path, '.fai')
downsample = opt$downsample
autosome_only = opt$autosome

width = 16 # inches
height = 7.2 # inches



counts_df = read_delim(counts_path, delim='\t', col_types = cols(chrom=col_character()))
reference_fai = read_delim(fai_path, delim='\t', col_names = c('name', 'length', 'offset', 'linebases', 'linewidth', 'qualoffset'))
reference_fai = reference_fai %>% mutate(cumpos = cumsum(length), lag1=lag(cumpos), lag1=replace_na(lag1, 0))


if(autosome_only==T){
  reference_fai = reference_fai %>% filter(str_detect(name, '^chr[0-9]+$|^[0-9]+$'))
  counts_df = counts_df %>% filter(str_detect(chrom, '^chr[0-9]+$|^[0-9]+$'))
}else{
  reference_fai = reference_fai %>% filter(str_detect(name, '^chr[0-9X]+$|^[0-9X]+$'))
  counts_df = counts_df %>% filter(str_detect(chrom, '^chr[0-9X]+$|^[0-9X]+$'))
  
}
total_genome_length = reference_fai %>% pull(length) %>% sum
reference_fai_dict = setNames(reference_fai$lag1, reference_fai$name)



chrom_to_numeric <- function(chromosome, position, reference_fai_dict){
  # change chromosomal position into a single numeric value for plotting
  return(reference_fai_dict[[chromosome]] + position)
}
chrom_to_numeric_vec = Vectorize(chrom_to_numeric, vectorize.args = c('chromosome', 'position'))


counts_df_numpos = counts_df %>% mutate(numeric_pos = chrom_to_numeric_vec(chrom, pos, reference_fai_dict))



plot_parent_vaf <- function(count_df, parent, reference_fai_dict, total_genome_length, downsample=1, plot_label='', print_chrom_label=F){
  parent_color = c('father'='#F9939B', 'mother'='#5882F9')
  color = parent_color[[parent]]
  pushViewport(viewport(x=0.05, y=0.2, width=0.9, height=0.77, just=c('left', 'bottom')))
  grid.rect()
  # grid.text(plot_label, x = 0, y=1.05, just=c('left', 'bottom'))
  for(i in 1:(length(reference_fai_dict))){
    fai_length = reference_fai_dict[i]
    if(i==length(reference_fai_dict)){
      next_length = total_genome_length
    }else{
      next_length = reference_fai_dict[i+1]
    }
    chrom_text_pos = (fai_length + next_length)/2
    
    grid.lines(x = c(fai_length/total_genome_length, fai_length/total_genome_length), y=c(0, 1))
    if(print_chrom_label==T){
      grid.text(label = gsub('chr', '', names(reference_fai_dict)[i]), x= chrom_text_pos/total_genome_length, y=-0.13, rot = 0, gp = gpar(fontsize = 14), just=c('center', 'center'))
      
    }
    
  }
  
  
  parent_vaf_column = paste0(parent, '_vaf')
  parent_depth_column = paste0(parent, '_depth')
  parent_alt_column = paste0(parent, '_alt')
  count_df = count_df %>% filter(get(parent_depth_column)-get(parent_alt_column) >1)
  
  count_df = count_df %>% filter(get(parent_vaf_column) <1)
  count_df = count_df %>% sample_frac(downsample)
  
  
  
  grid.yaxis(at=c(0, 0.5, 1), label=c(0, 0.5, 1), gp = gpar(fontsize = 14))
  # grid.text('VAF', rot = 90, x=-0.04, y=0.5, just=c('center', 'bottom'), gp = gpar(fontsize = 24))
  if(dim(count_df)[1]>0){
    grid.points(
      x=count_df$numeric_pos/total_genome_length, 
      y=count_df[[parent_vaf_column]], 
      gp=gpar(col=color, alpha=0.3), pch=16, size = unit(0.5, 'char')
    )
  }
  grid.lines(x=c(0, 1), y=c(0.5, 0.5), gp=gpar(col='black',lwd=2, lty=2))
  
  popViewport(1)
}

plot_parent_depth <- function(count_df, parent, reference_fai_dict, total_genome_length, downsample=1, print_chrom_label=F){
  pushViewport(viewport(x=0.05, y=0.2, width=0.9, height=0.72, just=c('left', 'bottom')))
  grid.rect()
  # grid.text('Depth', x = 0, y=1.05, just=c('left', 'bottom'), gp = gpar(fontsize = 18, fontface = "bold"))
  for(i in 1:(length(reference_fai_dict))){
    fai_length = reference_fai_dict[i]
    if(i==length(reference_fai_dict)){
      next_length = total_genome_length
    }else{
      next_length = reference_fai_dict[i+1]
    }
    chrom_text_pos = (fai_length + next_length)/2
    
    grid.lines(x = c(fai_length/total_genome_length, fai_length/total_genome_length), y=c(0, 1))
    if(print_chrom_label==T){
      grid.text(label =gsub('chr', '', names(reference_fai_dict)[i]), x= chrom_text_pos/total_genome_length, y=-0.13, rot = 0,  gp = gpar(fontsize = 14), just=c('center', 'center'))
    }
  }
  mean_depth = count_df %>% filter(!(str_detect(chrom, 'X|Y'))) %>% pull(depth) %>% mean(na.rm=T) %>% round
  
  count_df = count_df %>% sample_frac(downsample)
  
  parent_depth_column = paste0(parent, '_depth')
  
  if(dim(count_df)[1]>0){
    grid.points(
      x=count_df$numeric_pos/total_genome_length, 
      y=count_df[[parent_depth_column]]/(mean_depth*2), 
      gp=gpar(alpha=0.3), pch=16,  size = unit(0.5, 'char')
      )
  }
  grid.yaxis(at=seq(0, mean_depth*2, mean_depth)/(mean_depth*2), label=seq(0, mean_depth*2, mean_depth),  gp = gpar(fontsize = 14))
  # grid.text('Depth', rot = 90, x=-0.04, y=0.5, just=c('center', 'bottom'), gp = gpar(fontsize = 24))
  grid.lines(x=c(0, 1), y=c(0.5, 0.5), gp=gpar(col='gray',lwd=2, lty=2))
  
  popViewport(1)
}




# select homo-ref(mother)/homo-alt(father)
child_homoalt = counts_df_numpos %>% filter(vaf==1)

mother_contam_by_father = child_homoalt %>% filter(father_vaf==1)
father_contam_by_mother = child_homoalt %>% filter(mother_vaf==1)



grid.newpage()

if(opt$output_format=='png'){
  output_path = normalizePath(file.path(opt$output_dir, paste0(basename(counts_path), '.plot.png')))
  png(output_path, units='in', res=300, width=width, height=height)
}else if(opt$output_format=='jpg'){
  output_path = normalizePath(file.path(opt$output_dir, paste0(basename(counts_path), '.plot.jpg')))
  jpeg(output_path, units='in', res=300, width=width, height=height)
}else{
  output_path = normalizePath(file.path(opt$output_dir, paste0(basename(counts_path), '.plot.pdf')))
  pdf(output_path, width=width, height=height)
}

pushViewport(viewport(x=0.1, y=0, width=0.9, height=0.95, just=c('left', 'bottom')))
pushViewport(viewport(x=0.0,  y=0, width=1, height=1/6, just=c('left', 'bottom')))
plot_parent_vaf(father_contam_by_mother, 'father', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (mother) + het (father)', print_chrom_label=T)
popViewport(1)
pushViewport(viewport(x=0,  y=1/6, width=1, height=1/6, just=c('left', 'bottom')))
plot_parent_vaf(mother_contam_by_father, 'mother', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (mother) + het (father)', print_chrom_label=T)
popViewport(1)

pushViewport(viewport(x=0,  y=2/6, width=1, height=1/6, just=c('left', 'bottom')))
plot_parent_vaf(child_homoalt, 'father', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (mother) + homo-alt (father)')
popViewport(1)

pushViewport(viewport(x=0,  y=3/6, width=1, height=1/6, just=c('left', 'bottom')))
plot_parent_vaf(child_homoalt, 'mother', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (father) + homo-alt (mother)')
popViewport(1)


pushViewport(viewport(x=0,  y=4/6, width=1, height=0.2, just=c('left', 'bottom')))
plot_parent_depth(child_homoalt, 'father', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample)
popViewport(1)

pushViewport(viewport(x=0,  y=5/6, width=1, height=0.2, just=c('left', 'bottom')))
plot_parent_depth(child_homoalt, 'mother', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample)
popViewport(1)


popViewport(1)
pushViewport(viewport(x=0, y=0, width=0.1, height=1, just=c('left', 'bottom')))
grid.lines(x=c(0.99, 0.99), y=c(0.03, 2/6-0.03), gp=gpar(lwd=2)); grid.text('GroupE VAF', x=0.40, y=1/6, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 18)); grid.text('mother vaf=1', x=0.65, y=1/12, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 10));grid.text('father vaf=1', x=0.65, y=3/12, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 10));grid.text('father', x=0.85, y=1/12, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18)); grid.text('mother', x=0.85, y=3/12, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18))
grid.lines(x=c(0.99, 0.99), y=c(2/6+0.01, 4/6-0.03), gp=gpar(lwd=2)); grid.text('GroupD VAF', x=0.40, y=3/6, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 18)); grid.text('father', x=0.85, y=5/12, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18)); grid.text('mother', x=0.85, y=7/12, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18))
grid.lines(x=c(0.99, 0.99), y=c(4/6+0.01, 6/6-0.03), gp=gpar(lwd=2)); grid.text('Depth', x=0.40, y=5/6, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 18)); grid.text('father', x=0.85, y=9/12, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18)); grid.text('mother', x=0.85, y=11/12, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18))

dev.off()