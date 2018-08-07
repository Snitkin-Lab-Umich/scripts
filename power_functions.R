# FUNCTIONS TO CALCULATE POWER FOR BACTERIAL GWAS

# Functions included:
# min_snp_dists - to get minimum pairwise distances between pairs of isolates with different phenotypes
# power_locus - to calculate power by locus

# load libraries
library(ape)
library(pwr)

script_dir = dirname(sys.frame(1)$ofile)

source(paste0(script_dir,'/power_calculator_functions_farhat2014.R'))
source(paste0(script_dir,'/genomic_summary_functions.R'))

# Get minimum pairwise distances between pairs of isolates with different phenotypes
# Input:
# 1) aln - alignment (names are isolate_no)
# 2) aln_phen - alignment phenotype (same order as isolates in alignment)
# 3) phen_interest - character string of phenotype of interest in phen vector (ex. 'R' or 'yes')
# Output:
# 1) vector of minimum pairwise distances between each phenotype of interest and it's closest neighbor with a different phenotype
get_min_snp_dists = function(aln,aln_phen,phen_interest){
  # get minimum pairwise snp distances for phenotype pairs
  snpdist = dist.dna(aln,model = 'N',as.matrix = T)
  
  min_dists = 
    unlist(sapply(1:nrow(snpdist), function(x){
      if(aln_phen[x] == phen_interest){
        min(snpdist[x,][aln_phen != phen_interest])
      }else{
        NULL
      }
    }))
  
  return(min_dists)
}

# Function to calculate power of identifying associations in loci using pairs of isolates or chi-square test
# Input:
# 1) aln - alignment (to calculatepairwise snp distance)
# 2) tree - tree (to collapse isolates by monphyletic group)
# 3) phen - vector of phenotypes named by isolate_no
# 4) phen_interest - character string of phenotype of interest in phen vector (ex. 'R' or 'yes')
# 5) genome_length - numeric value corresponding to the genome length in bp
# 6) gene_lengths - vector of gene lengths
# 7) effect_size - list of effect sizes to test. Default: c(0.01,0.05,0.100,0.15,0.20,0.25,0.40,0.50,0.60)
# Note: Then names in the alignment and the names of the tip labels should be the isolate_no
# Output: list with:
# 1) pow - pairwise power
# 2) avg_snpdist_pairs - average snp distance for pairwise power
# 3) pow_subset - pairwise power after collapsing pure subtrees
# 4) avg_snpdist_pairs_subset - average snp distance for collapsed subtree power
# 5) pow_chisq - chi-square power results
# 6) es_chisq_80 - chi-square effect size for 80% power (for pairwise - right now, guess and check)
power_locus = function(aln,tree,phen,phen_interest,genome_length,gene_lengths,
                       effect_size = c(0.01,0.05,0.100,0.15,0.20,0.25,0.40,0.50,0.60)){
  
  # remove isolates with no phenotype associated with them
  no_phen = names(phen)[phen == '']
  
  aln = aln[!rownames(aln) %in% no_phen,]
  tree = drop.tip(tree,no_phen)
  
  aln_phen = sapply(rownames(aln), function(x) phen[names(phen) == x])
  names(aln_phen) = rownames(aln)
  
  tip_phen = sapply(tree$tip.label, function(x) phen[names(phen) == x])
  names(tip_phen) = tree$tip.label
  
  phen_names = names(table(phen))[names(table(phen)) != '']
  
  # get minimum pairwise snp distances for phenotype pairs
  snpdist = dist.dna(aln,model = 'N',as.matrix = T)
  
  min_dists = get_min_snp_dists(aln,aln_phen,phen_interest)
  
  # raw count (not correcting for clusters on tree)
  phen_count = sum(aln_phen == phen_interest)
  
  avg_snpdist_pairs = round(mean(min_dists),1)
  pow = powerCalculatorLocus(k = genome_length,s = avg_snpdist_pairs,f = effect_size,n = phen_count,gl = gene_lengths)
  
  # df_chisq = number of classes - 1
  pow_chisq = pwr.chisq.test(w = effect_size, N = length(tip_phen), df = 1)
  # effect size for power of 80
  es_chisq_80 = pwr.chisq.test(power = 0.8, w = NULL, N = length(tip_phen), df = 1)
  
  # collapse pure subtrees
  collapsed_tree = collapse_pure_subtrees(tree,tip_phen)
  
  aln_subset = aln[rownames(aln) %in% collapsed_tree$tip.label,]
  
  snpdist_subset = dist.dna(aln_subset,model = 'N',as.matrix = T)
  
  aln_phen_subset = sapply(rownames(aln_subset), function(x) phen[names(phen) == x])
  names(aln_phen_subset) = rownames(aln_subset)[rownames(aln) %in% collapsed_tree$tip.label]
  
  min_dists_subset = get_min_snp_dists(aln_subset,aln_phen_subset,phen_interest)
  
  phen_count_subset = sum(aln_phen_subset == phen_interest)
  
  avg_snpdist_pairs_subset = round(mean(min_dists_subset),1)
  pow_subset = powerCalculatorLocus(k = genome_length,s = avg_snpdist_pairs_subset,f = effect_size,n = phen_count_subset,gl = gene_lengths)
  
  return(list(pow=pow,
              avg_snpdist_pairs=avg_snpdist_pairs,
              pow_subset=pow_subset,
              avg_snpdist_pairs_subset=avg_snpdist_pairs_subset,
              pow_chisq=pow_chisq,
              es_chisq_80=es_chisq_80))
}

