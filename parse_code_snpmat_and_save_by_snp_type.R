# Katie Saund
# 2018-08-21

# This script will run quality control functions on Cdif specific snpmat. 

# Script takes in a raw code snp mat produced by Ali's variant calling pipeline. 

# SOURCE -----------------------------------------------------------------------
source("/nfs/esnitkin/Project_Cdiff/Analysis/Hanna_paper/2018-08-20_QC_snpmat/lib/2018-08-22_qc_snpmat_lib.R")

# PARSE ARGUMENTS --------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) #arguments from the PBS script

# code_snpmat should be a file: SNP_matrix_code.csv
code_snpmat <- read.csv(args[1], 
                        header = TRUE,
                        stringsAsFactors = FALSE,
                        sep = "\t",
                        quote = "", 
                        row.names = 1)

raw_phenotypes <- read.csv(args[2], 
                           row.names = 1, 
                           header = TRUE, 
                           sep = ",", 
                           stringsAsFactors = FALSE)
raw_phenotypes[raw_phenotypes == "ND"] <- NA

lookup <- read.table(args[3], 
                     header = TRUE, 
                     stringsAsFactors = FALSE)

tree <- read.tree(args[4])

# PARSE SNP MAT ----------------------------------------------------------------
simplified_code_snpmat <- simplify_snp_code(code_snpmat)
parsed_simple_code_snpmat <- parse_snps(simplified_code_snpmat)
parsed_simple_code_snpmat$snpmat <- fix_snpmat_isolate_names(code_snpmat, parsed_simple_code_snpmat$snpmat)
save(parsed_simple_code_snpmat, file = paste("../data/", format(Sys.time(), "%Y-%m-%d_"), "parsed_simple_code_snpmat.rda", sep = ""))

# CREATE MANHATTAN PLOTS OF SNP OCCURANCE --------------------------------------
create_manhattan_plot(parsed_simple_code_snpmat$pos, "all_snps")
create_manhattan_plot(parsed_simple_code_snpmat$pos[parsed_simple_code_snpmat$del_mut], "del_snps")
create_manhattan_plot(parsed_simple_code_snpmat$pos[parsed_simple_code_snpmat$ns_mut], "ns_snps")
create_manhattan_plot(parsed_simple_code_snpmat$pos[parsed_simple_code_snpmat$s_mut], "s_snps")

# PREPARE DATA FOR HEATMAPS ----------------------------------------------------
tree <- format_tree(tree, lookup)
tree <- drop.tip(tree, c(1:Ntip(tree))[!(tree$tip.label != "")])
ribo_info <- simplify_ribotype(raw_phenotypes, tree)
t_snpmat  <- fix_sample_names(parsed_simple_code_snpmat$snpmat, lookup, tree)
htmp_tree <- create_heatmap_compatible_tree(tree)

# gyrA & gyrB ------------------------------------------------------------------
# Want to find: gyrA C(245)T[T(82)I] and/or gyrB G(1276)A [D(426)N] 
# Should be almost exclusively in 027 isolates
gyr_mat <- subset_matrix_to_specific_genes(parsed_simple_code_snpmat, c("gyrA", "gyrB"), c(82, 426), t_snpmat)
save_heatmap(gyr_mat, htmp_tree, "fqR_mut_in_gyr", ribo_info$ribotype_color)

# TREHALOSE --------------------------------------------------------------------
trehalose_mat <- subset_matrix_to_specific_genes(parsed_simple_code_snpmat, "treR", c(171, 172), t_snpmat) 
save_heatmap(trehalose_mat, htmp_tree, "trehalose_mutations", ribo_info$ribotype_color)

# GET GENES WITH SPECIFIC TYPES OF MUTATIONS -----------------------------------
t_snpmat_with_gene_names <- rbind(t_snpmat, parsed_simple_code_snpmat$genes) 
rownames(t_snpmat_with_gene_names)[nrow(t_snpmat_with_gene_names)] <- "gene_names"
genes_with_stops       <- find_genes_with_specific_snp_types(parsed_simple_code_snpmat, parsed_simple_code_snpmat$stop,        t_snpmat, t_snpmat_with_gene_names, "genes_with_stops", "stop",                     htmp_tree, ribo_info$ribotype_color)
genes_with_ns_mut      <- find_genes_with_specific_snp_types(parsed_simple_code_snpmat, parsed_simple_code_snpmat$ns_mut,      t_snpmat, t_snpmat_with_gene_names, "genes_with_nonsense", "nonsense",              htmp_tree, ribo_info$ribotype_color)
genes_with_del_mut     <- find_genes_with_specific_snp_types(parsed_simple_code_snpmat, parsed_simple_code_snpmat$del_mut,     t_snpmat, t_snpmat_with_gene_names, "genes_with_deleterious", "deleterious",        htmp_tree, ribo_info$ribotype_color)
genes_with_snpeff_high <- find_genes_with_specific_snp_types(parsed_simple_code_snpmat, parsed_simple_code_snpmat$snpeff_high, t_snpmat, t_snpmat_with_gene_names, "genes_with_snpeff_high", "snpeff_high_impact", htmp_tree, ribo_info$ribotype_color)

# END OF SCRIPT ----------------------------------------------------------------