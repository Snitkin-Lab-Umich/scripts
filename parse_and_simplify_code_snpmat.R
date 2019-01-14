# Katie Saund
# 2019-01-14

# This script will (1) simplify a code snpmat then (2) parse it using the functions in the variant_parse_functions.R file. 

# Input: 
# 1. Path to code snp mat produced by Ali's variant calling pipeline. 
# 2. Path/filename of output snp mat. 
# Output: 
# 1. Simplified, parse binary snp mat. 

# SOURCE -----------------------------------------------------------------------
source("/nfs/esnitkin/bin_group/pipeline/Github/scripts/variant_parser_functions.R") # this is the version controlled script

# PARSE ARGUMENTS --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) # arguments from the PBS script

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
# parsed_simple_code_snpmat$snpmat <- fix_snpmat_isolate_names(code_snpmat, parsed_simple_code_snpmat$snpmat)
save(parsed_simple_code_snpmat, file = args[2])

# END OF SCRIPT ----------------------------------------------------------------