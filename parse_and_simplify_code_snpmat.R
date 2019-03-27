# Katie Saund
# 2019-01-14

# This script will (1) simplify a code snpmat then (2) parse it using the functions in the variant_parse_functions.R file. 

# Input: 
# 1. Path to code snp mat produced by Ali's variant calling pipeline. 

# Output: 
# 1. Simplified, parse binary snp mat as an rdata object. 

# SOURCE -----------------------------------------------------------------------
source("/nfs/esnitkin/bin_group/pipeline/Github/scripts/variant_parser_functions.R") # this is the version controlled script

# PARSE ARGUMENTS --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) # arguments from the PBS script

# Check arguments: 
snpmat_name <- "SNP_matrix_code.csv"
if (substr(args[1], nchar(args[1]) - nchar(snpmat_name) + 1, nchar(args[1])) != snpmat_name){
  stop("Name of input file should be SNP_matrix_code.csv")
}

# code_snpmat should be a file: SNP_matrix_code.csv
code_snpmat <- read.table(args[1], 
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t",
                        quote = "", 
                        row.names = 1)

# PARSE SNP MAT ----------------------------------------------------------------
simplified_code_snpmat <- simplify_snp_code(code_snpmat)
parsed_simple_code_snpmat <- parse_snps(simplified_code_snpmat)

# END OF SCRIPT ----------------------------------------------------------------