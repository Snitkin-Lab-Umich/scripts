# Katie Saund
# 2019-06-09
# Goal: this script reads in a snp_mat as created by Ali's pipeline and checks for previously reported bugs. 

library(tidyverse)
library(testthat)

snpmat <- read.table("/scratch/esnitkin_fluxod/apirani/Project_Cdiff/Analysis/Project_propensity_score_match/2019_05_14_PSM_variant_calling/core_temp_dir/test_data/SNP_matrix_allele_new.csv", 
                   header = TRUE,
                   stringsAsFactors = FALSE,
                   sep = "\t",
                   quote = "",
                   row.names = 1)

test_snp_mat_rownames <- function(smat){
  # Check that ERROR messages are not appearing
  expect_identical(grep("ERROR", row.names(smat)), integer(0))
  
  
  # Check that biopython "[342903:234239] bug is fixed
  expect_identical(grep("[[][0-9]", row.names(smat)), integer(0))
  
  # Non-coding snps should have two genes in the locus_tag (grep for the "-" dash)
  # Maybe be too strict because of issue with intragenic variants? TBD. 
  num_Non_Coding_rows <- sum((!is.na(str_match(row.names(smat), 'Non-Coding'))))
  num_intragenic <- sum((!is.na(str_match(row.names(smat), 'intragenic'))))
  num_Non_Coding_rows_with_two_genes <- length(row.names(smat)[(!is.na(str_match(row.names(smat), 'Non-Coding')))] %>% gsub(".*locus_tag=", "", .) %>% gsub(" Strand Information.*", "", .) %>% grep('-', .))
  expect_equal(num_Non_Coding_rows - num_intragenic, num_Non_Coding_rows_with_two_genes) 
  
  # What does this "WARNING_TRANSCRIPT_NO_START_CODON" mean? How should we handle it?
  expect_identical(grep("WARNING_TRANSCRIPT_NO_START_CODON", row.names(smat)), integer(0))
}

test_snp_mat_rownames(snpmat)


# testing on 2019-06-11 to fix a bug. 
temp <- row.names(snpmat)[!is.na(str_match(row.names(snpmat), 'Non-Coding'))]
foo <- temp[is.na(str_match(temp, 'intragenic'))]
summary(str_length(foo %>% gsub(".*locus_tag=", "", .) %>% gsub(" Strand Information.*", "", .)))
foo[(str_length(foo %>% gsub(".*locus_tag=", "", .) %>% gsub(" Strand Information.*", "", .))) == 11]

row.names(snpmat) %>% gsub(".*SNP at ", "", .) %>% gsub(" >.*", "", .) %>%  as.numeric() %>%  summary()
