# plot VAF, depth etc for triomix
# 2022.01.01 cjyoon draw whole genome wide vaf plot
# 2022.03.23 cjyoon draw with larger fonts for legibility
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
width = 16 # inches
height = 6 # inches
autosome_only = opt$autosome

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



plot_vaf <- function(count_df, parent, reference_fai_dict, total_genome_length, downsample=1, plot_label='', print_chrom_label=F){
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
  
  count_df = count_df %>% sample_frac(downsample)
  
  grid.yaxis(at=c(0, 0.5, 1), label=c(0, 0.5, 1), gp = gpar(fontsize = 14))
  # grid.text('VAF', rot = 90, x=-0.04, y=0.5, just=c('center', 'bottom'), gp = gpar(fontsize = 24))
  grid.points(
    x=count_df$numeric_pos/total_genome_length, 
    y=count_df$vaf, 
    gp=gpar(col=color, alpha=0.3), pch=16, size = unit(0.5, 'char')
  )
  grid.lines(x=c(0, 1), y=c(0.5, 0.5), gp=gpar(col='black',lwd=2, lty=2))
  
  popViewport(1)
}

plot_depth <- function(count_df, reference_fai_dict, total_genome_length, downsample=1, print_chrom_label=F){
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
  
  grid.points(
    x=count_df$numeric_pos/total_genome_length, 
    y=count_df$depth/(mean_depth*2), 
    gp=gpar(alpha=0.3), pch=16,  size = unit(0.5, 'char')
  )
  grid.yaxis(at=seq(0, mean_depth*2, mean_depth)/(mean_depth*2), label=seq(0, mean_depth*2, mean_depth),  gp = gpar(fontsize = 14))
  # grid.text('Depth', rot = 90, x=-0.04, y=0.5, just=c('center', 'bottom'), gp = gpar(fontsize = 24))
  grid.lines(x=c(0, 1), y=c(0.5, 0.5), gp=gpar(col='gray',lwd=2, lty=2))
  
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


pushViewport(viewport(x=0.1, y=0, width=0.9, height=1, just=c('left', 'bottom')))
pushViewport(viewport(x=0.0,  y=0, width=1, height=0.2, just=c('left', 'bottom')))
plot_vaf(homoref_het_father, 'father', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (mother) + het (father)', print_chrom_label=T)
popViewport(1)

pushViewport(viewport(x=0,  y=0.20, width=1, height=0.2, just=c('left', 'bottom')))
plot_vaf(homoref_het_mother, 'mother', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (father) + het (mother)')
popViewport(1)

pushViewport(viewport(x=0,  y=0.40, width=1, height=0.20, just=c('left', 'bottom')))
plot_vaf(homoref_homo_alt_father, 'father', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (mother) + homo-alt (father)')
popViewport(1)

pushViewport(viewport(x=0,  y=0.60, width=1, height=0.2, just=c('left', 'bottom')))
plot_vaf(homoref_homo_alt_mother, 'mother', reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample, plot_label='homo-ref (father) + homo-alt (mother)')
popViewport(1)


pushViewport(viewport(x=0,  y=0.80, width=1, height=0.2, just=c('left', 'bottom')))
plot_depth(bind_rows(homoref_het_father, homoref_het_mother, homoref_homo_alt_father, homoref_homo_alt_mother), reference_fai_dict, total_genome_length=total_genome_length, downsample=downsample)
popViewport(1)
popViewport(1)
pushViewport(viewport(x=0, y=0, width=0.1, height=1, just=c('left', 'bottom')))
grid.lines(x=c(0.99, 0.99), y=c(0.02, 0.39), gp=gpar(lwd=2)); grid.text('GroupB VAF', x=0.45, y=0.22, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 18)); grid.text('paternal', x=0.85, y=0.10, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18)); grid.text('maternal', x=0.85, y=0.30, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18))
grid.lines(x=c(0.99, 0.99), y=c(0.42, 0.79), gp=gpar(lwd=2)); grid.text('GroupA VAF', x=0.45, y=0.62, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 18)); grid.text('paternal', x=0.85, y=0.50, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18)); grid.text('maternal', x=0.85, y=0.70, just=c('center', 'center'), rot=90, gp=gpar(fontsize=18))
grid.lines(x=c(0.99, 0.99), y=c(0.82, 0.99), gp=gpar(lwd=2)); grid.text('depth', x=0.85, y=0.90, just=c('center', 'center'), rot=90,  gp = gpar(fontsize = 18))
popViewport(1)
dev.off()