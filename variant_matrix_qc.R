#' ---
#' params:
#'   mat: NULL
#'   title: Variant Matrix QC
#' ---
#' 
#' ---
#' title: `r params$title`
#' ---

#+ echo=F, warnings=T, message=T, error=T
# Input to script: SNP_matrix_code.csv OR Indel_matrix_code.csv

# Load libraries
# load packages
deps = c("pheatmap", "htmlTable", "filesstrings");
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE, repos = "http://cran.us.r-project.org", dependencies=TRUE);
  }
  suppressPackageStartupMessages(library(dep, verbose=FALSE, character.only=TRUE))
}

# Load functions
source('/nfs/esnitkin/bin_group/pipeline/Github/scripts/variant_parser_functions.R')

# params = list()
# params$mat = '/scratch/esnitkin_fluxod/apirani/Project_Penn_KPC/Analysis/Regional_KPC_transmission/2018_01_22_Penn_ST258_variant_calling/2018_12_19_10_12_11_core_results/data_matrix/matrices/Indel_matrix_code.csv'
# params$mat = '/scratch/esnitkin_fluxod/apirani/Project_Penn_KPC/Analysis/Regional_KPC_transmission/2018_01_22_Penn_ST258_variant_calling/2018_12_06_14_33_46_core_results/data_matrix/matrices/SNP_matrix_code.csv'

#' ### Input matrix:
#+ echo=F, warnings=T, message=T, error=F
params$mat

# Read in and parse variant matrix
if(grepl('.RData',params$mat)){
  load(params$mat)
  alt_mat = parsed
  }else{
  stop("Requires .RData output of variant parser as the input at this time")
  }
# }else if(grepl('SNP',params$mat)){
#   alt_mat = parse_snps(params$mat,save_rdata=T)
#   file.move(list.files('.','*RData'),'..')
# }else if(grepl('Indel',params$mat)){
#   alt_mat = parse_indels(params$mat,save_rdata=T)
#   file.move(list.files('.','*RData'),'..')
# }

print(length(parsed))
stop()
#' ### Notes
#' The input into gubbins is an alignment where all phage regions, repeat regions, manually masked regions, low FQ (heterozygous) variant positions, and low MQ (mapping quality) variant positions are masked. The gubbins output is then used to build a tree. The code variant matrices contain more information in that low FQ and low MQ variant positions are not masked, but the position in genomes where the variant was filtered are indicated with a -3 or -4, respectively. Below are plots of variant positions including and excluding low FQ and low MQ to determine whether or not to use these positions for downstream analyses such as GWAS.
#' 

#' ### What numbers mean: 
#+ echo=F, warnings=F, message=T, error=T

#' | Number | Description   |          
#' |:------:|:-------------:|
#' | -4     | filtered (low MQ) | 
#' | -3     | filtered (low FQ) | 
#' | -2     | filtered (phage/repeat region) | 
#' | -1     | unmapped | 
#' | 0      | reference allele |
#' | 1      | core variant  | 
#' | 2      | filtered (DP/QUAL) | 
#' | 3      | non-core variant |

#+ echo=F, warnings=F, message=T, error=T

#' ### Interpreting the heatmaps
#' Rows are genomes and columns are variants. The band at the top indicates the fraction of variants that are not filtered for some reason (the higher the better). The band on the left indicates the number of variants that are filtered in a given genome (the lower the better). Variants with a high fraction of filtered variants and genomes with a high number of filtered variants are suspect.
#' 
#'
#' ## Summary Statistics

#' Total number of variants:

#+ echo=F, warnings=F, message=T, error=T
nrow(alt_mat$code$mat)

#' Total number of non-phage variants:
#+ echo=F, warnings=F, message=T, error=T
sum(rowSums(alt_mat$code$mat == -2) == 0)
#+ echo=F, warnings=F, message=T, error=T
var_type_counts = (apply(as.data.frame(t(alt_mat$code$mat)), 2, function(x) table(factor(x,levels=c(-1,0,1,2,3,-2,-3,-4)))))
phage = var_type_counts[6,] != 0
lowfq = var_type_counts[7,] != 0 
lowmq = var_type_counts[8,] != 0 
unmapped = var_type_counts[1,] != 0
filtered = colSums(var_type_counts[c(4,7,8),]) != 0
true = var_type_counts[5,] != 0
only_true = var_type_counts[3,] != 0
var_count = ncol(var_type_counts)
nonphage_varcount = sum(rowSums(alt_mat$code$mat == -2) == 0)

#' Masked variant positions in alignment:
#+ echo=F, warnings=F, message=T, error=T
paste('Phage or repeat region (-2):', sum(phage), signif(sum(phage)/var_count,2))
paste('Low FQ (-3):', sum(lowfq),signif(sum(lowfq)/var_count,2))
paste('Low MQ (-4):', sum(lowmq),signif(sum(lowmq)/var_count,2))
#' Individual variant positions (not including phage):
#+ echo=F, warnings=F, message=T, error=T
paste('Unmapped, filtered, and true (-1, 2, 3):', sum(unmapped & filtered), signif(sum(unmapped & filtered)/nonphage_varcount,2))
paste('Only unmapped or true (-1, 3):', sum(unmapped & !filtered),signif(sum(unmapped & !filtered)/nonphage_varcount,2))
paste('Only filtered or true (2, 3):', sum(!unmapped & filtered),signif(sum(!unmapped & filtered)/nonphage_varcount,2))
paste('Only true (1):', sum(only_true),signif(sum(only_true)/nonphage_varcount,2))

# Heatmap function
heatmap_colors = c('#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788')
allele_heatmap = function(mat,col=heatmap_colors){
  if(dim(mat)[2] == 0){
    return('No variants fit this category.')
  }
  par(mfrow=c(1,1))
  pos_filt_frac = 1 - colSums(mat == 3 | mat == 1)/(colSums(mat == 2 | mat == -3 | mat == -4) + colSums(mat == 3 | mat == 1))
  names(pos_filt_frac) = colnames(mat)
  genome_filt_frac = rowSums(mat == 2 | mat == -3 | mat == -4)/ncol(mat)
  names(genome_filt_frac) = rownames(mat)
  pheatmap(mat,color=col,
           show_rownames = T, show_colnames = F,
           fontsize_row = 4,
           annotation_row = as.data.frame(genome_filt_frac), annotation_col = as.data.frame(pos_filt_frac),
           cluster_rows = T, cluster_cols = F)#,
           #cellwidth = 1/10, cellheight = 1/10,
           #width = 20, height = 12)
}

#' ## Phage Variants
#' 
#' Histogram of phage variants across the genome
#+ echo=F, warnings=F, message=T, error=T
hist(alt_mat$code$annots$pos[alt_mat$code$annots$phage],10000,main='',xlab='Phage variant positions')

#' ## All Variants
#' For all heatmaps: (1) before gubbins recombination filtering, (2) phage regions masked, (3) rows are genomes and columns are variants.
#' All variants (excludes phage).
#+ echo=F, warnings=F, message=T, error=T
mat = t(alt_mat$code$mat)
colnames(mat) = rownames(alt_mat$code$mat)
alt_mat$code$mat = NULL #to free up memory
print(table(colSums(mat == -2) == 0))
allele_heatmap(mat[,colSums(mat == -2) == 0],col=heatmap_colors[c(1,2,4:length(heatmap_colors))])

#' Only non-core variants.
#+ echo=F, warnings=F, message=T, error=T
allele_heatmap(mat[,filtered & !phage],col=heatmap_colors[c(1,2,4,5,7,8)])

#' ## Low FQ
#' Variants positions excluding those filtered because of low FQ.
#+ echo=F, warnings=T, message=T, error=T
# Mask variants with low FQ (-3) 
mat_maskFQ = mat
mat_maskFQ[,colSums(mat == -3) != 0] = -3
allele_heatmap(mat_maskFQ[,colSums(mat_maskFQ == -2) == 0 & colSums(mat_maskFQ == -3) == 0],col=heatmap_colors[c(1,2,4:length(heatmap_colors))])

#' Non-core variants that were filtered out because of low FQ. 
#+ echo=F, warnings=F, message=T, error=T
allele_heatmap(mat[,!phage & var_type_counts[7,] != 0],col=heatmap_colors[c(1,2,4:length(heatmap_colors))])

#' Non-core variants that were not filtered out because of low FQ. 
#+ echo=F, warnings=F, message=T, error=T
allele_heatmap(mat[,filtered & !phage & var_type_counts[7,] == 0],col=heatmap_colors[c(1,2,4:length(heatmap_colors))])

#' Histogram of low FQ positions across the genome
#+ echo=F, warnings=T, message=T, error=T
if(length(alt_mat$code$annots$pos[colSums(mat == -3) != 0]) != 0){
  hist(alt_mat$code$annots$pos[colSums(mat == -3) != 0],1000,main='',xlab='Low FQ positions across genome')
}

#' ## Low MQ
#' 
#' Variants positions excluding those filtered because of low FQ or low MQ.
#+ echo=F, warnings=T, message=T, error=T
# Mask variants with low MQ (-4) 
mat_maskMQ = mat_maskFQ
rm(mat_masFQ) #to reduce memory, since we don't need it any more
mat_maskMQ[,colSums(mat == -4) != 0] = -4
allele_heatmap(mat_maskMQ[,colSums(mat_maskMQ == -2) == 0 & colSums(mat_maskMQ == -3) == 0 & colSums(mat_maskMQ == -4) == 0],col=heatmap_colors[4:length(heatmap_colors)])

#' Non-core variants that were filtered out because of low MQ. 
#+ echo=F, warnings=F, message=T, error=T
allele_heatmap(mat[,!phage & var_type_counts[8,] != 0],col=heatmap_colors[c(1,2,4:length(heatmap_colors))])

#' Non-core variants that were not filtered out because of low FQ or low MQ. 
#+ echo=F, warnings=F, message=T, error=T
allele_heatmap(mat[,filtered & !phage & var_type_counts[7,] == 0 & var_type_counts[8,] == 0],col=heatmap_colors[c(4:(length(heatmap_colors)))])

#' Histogram of low MQ positions across the genome
#+ echo=F, warnings=T, message=T, error=T
if(length(alt_mat$code$annots$pos[colSums(mat == -4) != 0]) != 0){
hist(alt_mat$code$annots$pos[colSums(mat == -4) != 0],1000,main='',xlab='Low MQ positions across genome')
}

#' ## Variant count histograms.
#' 
#' For all these phage regions and variants positions filtered because of low FQ and low MQ have been removed.
#' 
#' Final variant counts across genome. 
#+ echo=F, warnings=F, message=T, error=T
hist(alt_mat$code$annots$pos,10000,col=rgb(1,0,0,1/4),border = rgb(1,0,0,1/4),main='',xlab='Variant positions')
if(length(alt_mat$code$annots$pos[colSums(mat_maskMQ == -2) == 0 & colSums(mat_maskMQ == -3) == 0 & colSums(mat_maskMQ == -4) == 0]) != 0){
hist(alt_mat$code$annots$pos[colSums(mat_maskMQ == -2) == 0 & colSums(mat_maskMQ == -3) == 0 & colSums(mat_maskMQ == -4) == 0],10000,add=T,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/4))
legend('topright',c('Pre-masking','Post-masking\n(input to gubbins)'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)
}
rm(alt_mat) #no longer need; to reduce memory

#' Number of genomes with each variant. 
#+ echo=F, warnings=T, message=T, error=T
hist(colSums(mat == 3 | mat == 1),1000,main='',xlab='Number of genomes with variant',col=rgb(1,0,0,1/4),border = rgb(1,0,0,1/4))
hist(colSums(mat_maskMQ == 3 | mat_maskMQ == 1),1000,main='',xlab='Number of genomes with variant',add=T,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/4))
legend('topright',c('Pre-masking','Post-masking\n(input to gubbins)'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)


#' Number of variants per isolate.
#+ echo=F, warnings=T, message=T, error=T
barplot(sort(rowSums(mat == 3 | mat == 1),decreasing = T),1000,main='',xaxt='n',xlab='Genomes',ylab='Numer of variants',col=rgb(1,0,0,1/4),border = rgb(1,0,0,1/4))
barplot(sort(rowSums(mat_maskMQ == 3 | mat_maskMQ == 1),decreasing = T),1000,main='',xaxt='n',xlab='Genomes',ylab='Numer of variants',add=T,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/4))
legend('topright',c('Pre-masking','Post-masking\n(input to gubbins)'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)

#' Number of uniquely filtered positions in genome
#+ echo=F, warnings=T, message=T, error=T
binmat = ((mat == 2) | (mat == -3) | (mat == -4))
uniq_filt_ct = (rowSums(binmat[,which(colSums(binmat) == 1)]))
hist(uniq_filt_ct,1000,main='',xlab='Number of uniquely filtered positions in genome',col=rgb(1,0,0,1/4),border = rgb(1,0,0,1/4))
rm(mat) # no longer need, to free up memory
binmat = ((mat_maskMQ == 2) | (mat_maskMQ == -3) | (mat_maskMQ == -4))
uniq_filt_ct = (rowSums(binmat[,which(colSums(binmat) == 1)]))
hist(uniq_filt_ct,1000,main='',xlab='Number of uniquely filtered positions in genome',add=T,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/4))
legend('topright',c('Pre-masking','Post-masking\n(input to gubbins)'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)
