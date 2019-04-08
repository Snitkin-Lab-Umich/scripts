# Functions to parse panISa output
# panISa on GitHub: https://github.com/bvalot/panISa

# Load libraries
library(genbankr)

# Function to parse is results into a matrix. Takes the output of the panisa isfinder script as input.
# Input: ismat: panISa isfinder script output text file
# Output: matrix where rows are IS element named by start_end position and columns are sample. 1 indicates presence, 0 indicates absence
make_ismat = function(isfinder){
  if(is.character(isfinder)){
    isfinder = read.delim(isfinder)
  }
  isfinder = isfinder[!grepl('Sample',isfinder$Sample),]
  uniq_samp = unique(isfinder$Sample)
  pos = paste0(isfinder$Start_Position,'_',isfinder$Stop_Position)
  uniq_pos = unique(pos)
  ismat = matrix(0,nrow=length(uniq_pos),ncol=length(uniq_samp))
  rownames(ismat) = uniq_pos
  colnames(ismat) = uniq_samp
  for(i in 1:nrow(isfinder)){
    ismat[pos[i],isfinder$Sample[i]] = 1
  }
  colnames(ismat) = gsub('_panISa','',colnames(ismat))
  return(ismat)
}

# Function to annotate row names of the IS matrix (made with parse_is) to include locus_tag information
# Input: ismat: insertion matrix from parse_is
#        gb: genbank file of a reference genome
# Output: matrix where row names contain locus_tag information. 
#         & indicates the IS element spans more than 1 gene
#         | indicates that the IS element is intergenic
annotate_is = function(ismat,gb){
  if(is.character(gb)){
    gb = readGenBank(gb)
  }
  gene_info = genes(gb)
  locus_tag = gene_info$locus_tag
  ranges = ranges(gene_info)
  gene_strand = as.character(strand(gene_info))
  is_locus_tags = sapply(rownames(ismat), function(x){
    is_pos = strsplit(x,'_')[[1]]
    is_start = as.numeric(is_pos[1])
    is_end = as.numeric(is_pos[2])
    is_lt_start = locus_tag[is_start >= start(ranges) & is_start <= end(ranges)]
    is_lt_end = locus_tag[is_end >= start(ranges) &  is_end <= end(ranges)]
    ig = sum(length(is_lt_start),length(is_lt_end)) == 0
    if(ig){
      ind = which(is_start >= start(ranges) & is_start <= lead(end(ranges)) | is_end >= start(ranges) & is_end <= lead(end(ranges)))
      is_lt_igs = locus_tag[c(ind,ind+1)]
      is_lt_ig = paste(is_lt_igs,collapse = '-')
    }else{
      is_lt_igs = NULL
      is_lt_ig = NULL
    }
    strand = gene_strand[locus_tag %in% c(is_lt_start,is_lt_igs,is_lt_end)]
    if(ig){
      paste0(is_lt_ig,'|',paste(gene_strand[c(ind,ind+1)],collapse='.')) # intergenic
    }else{
      paste0(paste(unique(c(is_lt_start,is_lt_end)),collapse='&'),'|',paste(strand,collapse='&')) # & means in more than 1 gene
    }
  })
  rownames(ismat) = paste0(rownames(ismat),'|',is_locus_tags)
  return(ismat)
}

parse_is = function(ismat){
  # SEPARATE ANNOTATIONS
  all_info = strsplit(rownames(ismat),'\\|')
  
  # GET START AND END
  start_end = sapply(all_info, function(x) x[[1]])
  start = gsub('_.*','',start_end)
  end = gsub('.*_','',start_end)
  
  # GET LOCUS TAGS
  locus_tag = sapply(all_info, function(x) x[[2]])
  # INTERGENIC
  locus_tag_ig_gene1 = sapply(strsplit(locus_tag,'-|&'),function(x) x[1])
  locus_tag_ig_gene2 = sapply(strsplit(locus_tag,'-'),function(x) x[2])
  # SPANS TWO GENES
  locus_tag_gene1 = sapply(strsplit(locus_tag,'-|&'),function(x) x[1])
  locus_tag_gene2 = sapply(strsplit(locus_tag,'&'),function(x) x[2])
  
  # INTERGENIC BOOLEAN
  intergenic = !is.na(locus_tag_ig_gene2)
  # SPANS TWO GENES BOOLEAN
  two_genes = !is.na(locus_tag_gene2)
  
  # GET STRAND 
  strand = sapply(all_info, function(x) x[[3]])
  # INTERGENIC
  strand_ig_gene1 = sapply(strsplit(strand,'\\.|&'),function(x) x[1])
  strand_ig_gene2 = sapply(strsplit(strand,'\\.'),function(x) x[2])
  # SPANS TWO GENES
  strand_gene1 = sapply(strsplit(strand,'\\.|&'),function(x) x[1])
  strand_gene2 = sapply(strsplit(strand,'&'),function(x) x[2])
  
  # INTERGENIC UPSTREAM
  upstream_ig_gene1 = intergenic & strand_ig_gene1 == '-'
  upstream_ig_gene2 = intergenic & strand_ig_gene2 == '+'
  
  high_impact = !intergenic | upstream_ig_gene1 | upstream_ig_gene2
  
  parsed = list(mat = ismat,
                start = start,
                end = end,
                locus_tag = locus_tag,
                locus_tag_ig_gene1 = locus_tag_ig_gene1,
                locus_tag_ig_gene2 = locus_tag_ig_gene2,
                locus_tag_gene1 = locus_tag_gene1,
                locus_tag_gene2 = locus_tag_gene2,
                intergenic = intergenic,
                two_genes = two_genes,
                strand = strand,
                strand_ig_gene1 = strand_ig_gene1,
                strand_ig_gene2 = strand_ig_gene2,
                strand_gene1 = strand_gene1,
                strand_gene2 = strand_gene2,
                upstream_ig_gene1 = upstream_ig_gene1,
                usptream_ig_gene2 = upstream_ig_gene2,
                high_impact = high_impact
                )
  
  return(parsed)
}

