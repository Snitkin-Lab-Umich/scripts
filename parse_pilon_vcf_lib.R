# Katie Saund  
# 2018-07-04  
#
# This library has functions to parse pilon multi-vcf output.  
# The input file for this parser is pilon vcf results that have been sorted, 
# indexed, and merged into one vcf. Use a version of this PBS script:
# /nfs/esnitkin/Project_Cdiff/Analysis/lib/2018-07-04_format_pilon_variants/2018-07-04_merge_pilon_vcfs.pbs
# to compile such a file from your pilon results. Remember to update copy 
# of the PBS with your username and flux allocation. 
#
# An example of using these functions is shown in: 
# /nfs/esnitkin/Project_Cdiff/Analysis/lib/2018-07-04_format_pilon_variants/2018-07-04_example_of_parsing_pilon_vcf_results.R
# Output of this parse should be 2 sets of data: large variants and snps. 
# a) pilon_large_variants
#   i) $allele_matrix 
#     (row names = cat(pos, svtype, svlen)) -> indicates which alt allele (0/1/2/3/etc...)
#   ii) $svtype
#     char = "INS" or "DEL" 
#   iii) $svlen
#     num
#   iv) $pos
#     num
#   v) $ref
#     char
#   vi) $alt 
#     list of chars
#   vii) $info
#     list of chars
#   viii) $imprecise
#     char == "IMPRECISE" or "PRECISE"
# 
# b) pilon_snps
#   i) $allele_matrix 
#     (row names = cat(pos, svtype, svlen)) -> indicates which alt allele (0/1/2/3/etc...)
#   ii) $pos
#     num
#   iii) $ref
#     char
#   iv) $alt 
#     list of chars
#   v) $info
#     list of chars

# PACKAGES ---------------------------------------------------------------------
library(data.table)

# LIBRARY ----------------------------------------------------------------------
report_data_summary <- function(vcf){
  print(paste(nrow(vcf), "variants identified", sep = " "))
  print(paste(ncol(vcf) - 9, "samples in the vcf file", sep = " "))
  print(paste(100 * sum(vcf$FILTER == "PASS") / nrow(vcf), "% of variants PASS the default filter check"), sep = "")
  print(paste(unique(vcf$'#CHROM'), "is the reference genome"))
} # end report_data_summary()

simplify_gt_calls <- function(vcf){
  # get all sample gt columns:
  vcf_colnames = colnames(vcf)
  column_index <- c(1:ncol(vcf))
  sample_column_index = column_index[!(colnames(vcf) %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))]
  
  # Simplifies the genotype (gt) calls in the VCF.
  vcf_info = vcf[,1:9]
  vcf_samples = apply(vcf[,sample_column_index], 2, function(x){
    n = gsub('/.*','',x)
    n = gsub('\\.','0',n)
    as.numeric(n)
  })
  vcf = as.data.frame(cbind(vcf_info,vcf_samples))
  colnames(vcf) = vcf_colnames
  
  if(length(grep("/", unique(unlist(vcf[ , sample_column_index])))) > 0){
    stop("Genotype alleles include ambiguous calls. Repeat vcf merge to filter out ambiguous results.")
  } 
  return(vcf)
} # end simplify_gt_calls()

get_variant_lengths <- function(vcf){
  # Determine length of allele 
  var_len <- rep(list(), nrow(vcf))
  for (i in 1:nrow(vcf)){
    temp <- strsplit(vcf$ALT[i], split =  ",")
    temp_len <- rep(0, length(temp[[1]]))
    for (j in 1:length(temp[[1]])){
      temp_len[j] <- nchar(temp[[1]][j])
    }
    var_len[[i]] <- temp_len
  }
  return(var_len)
} # end get_variant_lengths()

get_ref_length <- function(vcf){
  # Determine length of reference
  ref_len <- rep(0, nrow(vcf))
  for (i in 1:nrow(vcf)){
    ref_len[i] <- nchar(vcf$REF[i])
  }
  return(ref_len)
} # end get_ref_length()

subset_to_snps <- function(vcf){
  # Get just SNPs from the VCF (no long indels)
  var_length <- get_variant_lengths(vcf)
  ref_length <- get_ref_length(vcf)
  snp_index <- matrix(0, nrow = 1, ncol = 2)
  colnames(snp_index) <- c("row", "alt")
  for (i in 1:length(var_length)){
    for (j in 1:length(var_length[[i]])){
      if(ref_length[i] == 1 & var_length[[i]][j] == 1){
        snp_index <- rbind(snp_index, c(i , j))
      }
    }
  }
  snp_index <- snp_index[2:nrow(snp_index), ]
  all_snp <- vcf[unique(snp_index[ , 1]), ]
  return(all_snp)
} # end subset_to_snps()

subset_to_indels_length_n <- function(vcf, indel_length){
  # Get just the indels of at length of at least indel_length. 
  var_length <- get_variant_lengths(vcf)
  ref_length <- get_ref_length(vcf)
  indel_index <- matrix(0, nrow = 1, ncol = 2)
  colnames(indel_index) <- c("row", "alt")
  for (i in 1:length(var_length)){
    for (j in 1:length(var_length[[i]])){
      if(var_length[[i]][j] >= indel_length | ref_length[i] >= indel_length){
        indel_index <- rbind(indel_index, c(i , j))
      }
    }
  }
  indel_index <- indel_index[2:nrow(indel_index), ]
  all_indel <- vcf[unique(indel_index[ , 1]), ]
  return(all_indel)
} # end subset_to_indels_length_n()

convert_vcf_to_allele_matrix <- function(vcf){
  # Subset VCF to just the allele matrix. 
  allele_matrix <- vcf
  for (i in 1:nrow(allele_matrix)){
    row.names(allele_matrix)[i] <- paste(vcf$POS[i], ";", 
                                         vcf$REF[i], ";",
                                         vcf$ALT[i], ";", 
                                         sep = "")
  }
  column_index <- c(1:ncol(vcf))
  sample_column_index <- column_index[!(colnames(vcf) %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))]
  allele_matrix <- allele_matrix[ , sample_column_index]
  return(allele_matrix)
} # end convert_vcf_to_allele_matrix()

pilon_parse_snps <- function(vcf){
  # Convert VCF into a SNP matrix with associated information. 
  # b) pilon_snps
  #   i) $allele_matrix 
  #         (row names = cat(pos, svtype, svlen)) -> indicates which alt allele (0/1/2/3/etc...)
  #   ii) $pos
  #         num
  #   iii) $ref
  #         char
  #   iv) $alt 
  #         list of chars
  #   v) $info
  #         list of chars
  
  all_snp_results <- subset_to_snps(vcf)
  allele_matrix <- convert_vcf_to_allele_matrix(all_snp_results)
  pos <- all_snp_results$POS
  ref <- all_snp_results$REF
  alt <- rep(list(), nrow(all_snp_results))
  for (i in 1:nrow(all_snp_results)){
    alt[[i]] <- strsplit(all_snp_results$ALT[i], split =  ",")
  }
  info <- all_snp_results$INFO
  stopifnot(table(length(alt)) == 1)
  stopifnot(table(length(ref)) == 1)
  return(list("allele_matrix" = allele_matrix, 
              "pos" = pos,
              "ref" = ref, 
              "alt" = alt,
              "info" = info))
} # end pilon_parse_snps()

pilon_parse_large_variants <- function(vcf, min_indel_size){
  # Convert VCF into an indel matrix with associated information. 
  # Will only include indels of at size equal to or greater than indel_size 
  # a) pilon_large_variants
  #   i) $allele_matrix 
  #         (row names = cat(pos, svtype, svlen)) -> indicates which alt allele (0/1/2/3/etc...)
  #   ii) $svtype
  #         char = "INS" or "DEL" 
  #   iii) $svlen
  #         num
  #   iv) $pos
  #         num
  #   v) $ref
  #         char
  #   vi) $alt 
  #         list of chars
  #   vii) $info
  #         list of chars
  #   viii) $imprecise
  #         char == "IMPRECISE" or "PRECISE"
  
  long_indels <- subset_to_indels_length_n(vcf, min_indel_size)
  allele_matrix <- convert_vcf_to_allele_matrix(long_indels)
  svtype <- svlen <- NULL
  temp_info <- strsplit(long_indels$INFO,";")
  for (i in 1:nrow(long_indels)){
    for (j in 1:length(temp_info[[i]])){
      if (length(grep("SVLEN", temp_info[[i]][j]))){
        svlen[i] <- as.numeric(gsub("SVLEN=", "", temp_info[[i]][j]))
      } 
      if (length(grep("SVTYPE", temp_info[[i]][j]))){
        svtype[i] <- as.character(gsub("SVTYPE=", "", temp_info[[i]][j]))
      } 
    }
  }
  pos <- long_indels$POS
  ref <- long_indels$REF
  alt <- rep(list(), nrow(long_indels))
  for (i in 1:nrow(long_indels)){
    alt[[i]] <- strsplit(long_indels$ALT[i], split =  ",")
  }
  info <- long_indels$INFO
  imprecise <- rep(FALSE, length(info))
  imprecise[grep("IMPRECISE", info)] <- TRUE
  
  return(list("allele_matrix" = allele_matrix, 
              "svtype"= svtype,
              "svlen" = svlen,
              "pos" = pos,
              "ref" = ref, 
              "alt" = alt,
              "info" = info,
              "imprecise" = imprecise))
} # end pilon_parse_large_variants()
