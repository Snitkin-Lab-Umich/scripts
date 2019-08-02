# Functions to rereference variants based on a tree with an outgroup

# Load libraries
library(ape)
library(tidyverse)

# Load functions
source('/nfs/esnitkin/bin_group/pipeline/Github/scripts/variant_parser_functions.R')

# Combine snp and indel matrices
# Input:
# (1) snpmat: path to snp matrix file or rdata file, or parsed output itself
# (2) indelmat: path to indel matrix file or rdata file, or parsed output itself
# Output: 
# Combined snp and indel matrix of 0s and 1s
comb_snp_indel_mat = function(snpmat,indelmat){
  # Load in snpmat
  if(is.character(snpmat)){
    #if(grep('.RData',snpmat)){
    if(!is.na(str_match(snpmat,'.RData'))){
      load(snpmat)
      snps_parsed = parsed
    } else{
      snps_parsed = parse_snps(snpmat)
    }
  }else{
    snps_parsed = parse_snps(snpmat)
  }
  # Load in indelmat
  if(is.character(indelmat)){
    #if(grep('.RData',indelmat)){
    if(!is.na(str_match(indelmat,'.RData'))){
      load(indelmat)
      indels_parsed = parsed
    }else{
      indels_parsed = parse_snps(indelmat)
    }
  }else{
    indels_parsed = parse_snps(indelmat)
  }
  # Remove loaded parsed file (already stored in another variable)
  rm(parsed)
  # Combine snp and indel matrices
  mat = rbind(snps_parsed$mat,indels_parsed$mat)
  # Change matrix to zeros and ones
  mat[mat == 1 | mat == 3] = 1
  mat[mat != 1] = 0 
  # Fix column names
  colnames(mat) = gsub('_R1.fastq.gz','',colnames(mat))
  colnames(mat) = gsub('S.*_L*1.fastq.gz','',colnames(mat))
  colnames(mat) = gsub('_R1_001.fastq.gz','',colnames(mat))
  colnames(mat) = gsub('_1.fastq.gz','',colnames(mat))
  # Keep only rows with variants (but not all variants )
  mat = mat[rowSums(mat) > 0 & rowSums(mat) != ncol(mat),]
  return(mat)
}

# Root tree based on outgroup
# Input:
# (1) tree: file to tree (with outgroup) to be used for rereferencing, or tree itself
# (2) outgroup: name of outgroup in tree (default: 'outgroup')
# Output:
# rooted tree without outgroup
root_tree_og = function(tree,outgroup='outgroup'){
  if(is.character(tree)){
    # Load in tree
    tree = read.tree(tree)
  }
  # Root tree on outgroup
  tree = root(tree,outgroup)
  tree = drop.tip(tree,outgroup)
  return(tree)
}

# Get ancestral state of root
# Input:
# (1) tree - rooted tree to perform ancestral reconstruction on
# (2) tip_states - vector of tip states in same order as tree tip labels (if named, will check order is correct)
get_root_anc_state = function(tree,tip_states){
  if(!is.null(names(tip_states))){
    # Check if tree and tip_states match
    if(sum(tree$tip.label != names(tip_states)) > 0){
      stop('Tree tip labels and tip_state names do not match.')
    }
  }
  # If all variants, prob ancestral variant is 1 = 1
  if(sum(tip_states) == length(tip_states)){
    return(c(0,1))
  }
  # If none variants, prob ancestral variant is 0 = 1
  if(sum(tip_states) == 0){
    return(c(1,0))
  }
  cts = table(tip_states)
  # If only one 0 or one 1, ancestral state is state that occurs more than once
  if(length(cts[cts>1])==1){
    ar = c(0,0)
    ar[which.max(cts)] = 1
    return(ar)
  }
  # If ancestral state not obvious, perform ancestral reconstruction
  ar = ace(x = tip_states,phy = tree,type = 'discrete')
  # Probability that ancestral state of root is 0 or 1
  return(ar$lik.anc[1,])
}


# Rereference variants based on tree with an outgroup
# Input:
# (1) tree: rooted tree
# (2) variant matrix of 0s and 1s (rows are variants, columns are samples) (columns in same order as the tree tip labels)
# Output:
# (1) Returns rereferenced matrix
# (2) Saves .RData file ([date]_reref_mat.RData)
reref_vars = function(tree,mat,mattype='varmat',root_prob_hist=T){
  
  if(sum(!(tree$tip.label %in% colnames(mat))) > 0){
    stop('Some tree tip labels not in varinat matrix.')
  }
  
  # Keep only columns in tree (and order correctly)
  mat = mat[,tree$tip.label]
  
  if(sum(tree$edge.length == 0) > 0){
    warning('All zero branch lengths changed to small non-zero number to be able to perform ancestral reconstruction.')
    # Change any edge lengths that are zero to a very small number (so ancestral reconstruction doesn't break)
    tree$edge.length[tree$edge.length == 0] = min(tree$edge.length[tree$edge.length > 0])/1000
  }
  
  
  
  # Get ancestral state of root; 1st column = var absent (0), 2nd column = var present (1)
  ar_all = apply(mat,1,function(tip_states){
    get_root_anc_state(tree,tip_states)
  })
  
  if(root_prob_hist){
    dir.create('figures',showWarnings = F)
    file=paste0('figures/',format(Sys.time(),'%Y-%m-%d'),'_anc_state_prob_hist_',mattype,'.pdf')
    pdf(file)
    hist(ar_all[1,],100,xlab='Probability that ancestral state is 0 (before re-referencing)',main='')
    dev.off()
  }
  
  # Get most probable ancestral state
  ar_bin = apply(ar_all,2,which.max)
  ar_bin = ar_bin - 1
  
  # Re-reference matrix based on root ancestral state
  reref_mat = t(sapply(1:length(ar_bin), function(x){
    if(ar_bin[x] == 0){
      # Don't change if 0 is already the ancestral state
      mat[x,]
    }else{
      # Switch 0s and 1s if 1 is the ancestral state
      as.numeric(!mat[x,])
    }
  }))
  
  rownames(reref_mat) = rownames(mat)
  
  # Save re-referenced matrix
  save(reref_mat,file=paste0(format(Sys.time(),'%Y-%m-%d'),'_reref_',mattype,'.RData'))
  
  # Return re-referenced matrix
  return(reref_mat)
}

