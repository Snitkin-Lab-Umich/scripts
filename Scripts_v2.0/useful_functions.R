#Useful libraries
library(ape)
library(phytools)
library(seqinr)
library(RColorBrewer)
library(gplots)
library(phangorn)
library(stringr)
library(ggtree)
library(BiocInstaller)
library(ggbio)
library(Biostrings)



plotphylo = function(tree, clustervar, type="fan", cex=0.75, pch=16, col = c(colors3, colors200), main="", legend = F, columnlength = 10, sub="", legend_title="",
                     edge.width = 1){
  #clustervar is a named integer vector where the names are the tip labels and the ints are the cluster number
  if(legend == FALSE){
    plot(tree, type=type, show.tip.label=F, main=main, sub = sub, edge.width = edge.width)
    tips = which(tree$tip.label %in% names(clustervar))
    clustervar = clustervar[names(clustervar) %in% tree$tip.label]
    tiplabels(tip = tips, cex = cex, pch=pch, 
              col = col[clustervar[order(match(names(clustervar), tree$tip.label))]])
  }
  if (legend== TRUE){
    par(oma = c(1,1,1,10)) #makes plot with space on the right (large margin)
    #If you want space on the bottom make it oma = c(10,1,1,1) - goes (bottom, left, top, right)
    #You can also play around with the 10 and make it bigger or smaller
    plotphylo(tree, clustervar, type, cex, pch, col, main, legend = FALSE, columnlength, sub, legend_title, edge.width = edge.width) #plot main figure
    par(fig = c(0, 1, 0, 1), oma = c(1, 1, 1, 1), mar = c(0, 0, 0, 0), new = TRUE) #make new margins the whole space
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n") #plot blank figure on whole space
    #add legend or whatever you want to the right. 
    legend(x="right", horiz = F, legend=(sort(unique(clustervar))), 
           bty="n", col=col, pch=16, cex = 1, x.intersp = 1, y.intersp = 0.75,
           ncol = ceiling(length(unique(clustervar))/columnlength), title=legend_title)
  }
}

plot_clusters = function(cl_list, tree, titles=NULL, cex = 0.75, type = "phylogram"){
  #n = length(tree$tip.label)
  sapply(1:length(cl_list), function(i){
    x = cl_list[[i]]
    #clustervar = c(rep(1,length(x)), rep(2, n-length(x)))
    #not_in_cluster = tree$tip.label[which(!(tree$tip.label %in% x))]
    #names(clustervar) = c(x,not_in_cluster)
    #plotphylo(tree, clustervar, col = c("red", "grey"), type = type, main = titles[i], cex = cex)
    tips = match(x, tree$tip.label)
    plot(tree, show.tip.label = F, type = type, main = titles[i], edge.width = 0.5)
    tiplabels(tip = tips, cex = cex, pch = 16, col = "red")
  })
}


####NJ CLUSTERS#####
all_NJclusters <- function(vcf, min = 3){
  # calculates how many unique variants are in a cluster
  # inputs: 
  # vcf: matrix with columns of genomes and rows of variant positions
  # min: minimum acceptable number of shared variants for a cluster to be returned
  
  #clean vcf if not already 
  if("ALT" %in% colnames(vcf)){
    #get rid of variants whose reference or alternative allele is N
    vcf = vcf[!(grepl("N", vcf$ALT)),] 
    vcf = vcf[!(grepl("N", vcf$REF)),]
    reps = sapply(vcf$POS, function(x){
      return(length(which(vcf$POS == x)))
    })
    vcf = vcf[which(reps==1),]
    rownames(vcf) = vcf$POS #make variant positions row names
    vcf = vcf[,10:ncol(vcf)] #get rid of explanatory columns, just keep genome names
  }
  basemat = base_matrix(vcf)
  isolates = colnames(vcf) 
  clusters = vector(mode="list", length = length(isolates))
  
  for(i in 1:length(isolates)){
    write(i, file="/nfs/esnitkin/Project_Cdiff/Analysis/targeted_sequencing_paper/2018-05-16_NJ_cluster_picker/2018-07-18_RM_variants/isolates", append = T) 
    cls = NJcluster(vcf, min, basemat, start = isolates[i], clusters)
    if(length(cls)>0){
      n = first_null(clusters)
      clusters[n:(n+length(cls)-1)] <- cls
    }
  }
  
  #find and return unique clusters
  l = clusters[-which(sapply(clusters, is.null))]
  l = sapply(l, function(x){
    x = sort(x)
  })
  l = unique(l)
  return(l)
  
}


first_null = function(list){
  for(i in 1:length(list)){
    if(is.null(list[[i]])){
      return(i)
    }
  }
  return(length(list)+1)
}

NJcluster <- function(vcf, min = 3, basemat, start = NULL, set_clusters){
  # finds clusters that have unique defining variants that all include a specific genome
  # inputs: 
  # vcf: matrix with columns of genomes and rows of variant positions
  # min: minimum acceptable number of shared variants for a cluster to be returned
  # basemat: matrix with the same number of rows as vcf and 4 columns which are the number of A's, T's, G's & C's in that position
  # start: genome to start clustering at- if NULL, clustering starts with 2 most related genomes in vcf
  # set_clusters: clusters already defined
  clusters = vector(mode="list", length = (ncol(vcf) - 1)) #initialize list of possible clusters
  num_unique_vars = numeric(ncol(vcf) - 1) #initialize vector of number of unique variants for each cluster
  set_clusters = sapply(set_clusters, function(x){
    x = sort(x)
  })
  if(is.null(start)){ #if start is not specified
    cluster = two_closest(vcf)  #pick 2 isolates that are closest
  }
  else{ #if start is specified
    add = find_closest(vcf, start) #find closest other genome to specified starting genome
    cluster = c(start,add) #make those 2 genomes the working cluster 
  }
  if(in_list(sort(cluster),set_clusters)){ #if cluster is already picked, break out of loop
    return(list())
  }
  clusters[[1]] <- cluster #store cluster
  
  num_unique_vars[1] = num_unique_variants(vcf, cluster, basemat) #calculate number of unique variants, store
  #print(num_unique_vars)
  for(i in 2:(length(num_unique_vars))){ #keep adding one genome at a time until all are in the cluster
    reduced_vcf = reduce_vcf(vcf, cluster) #reduce vcf to just variants that are the same in cluster
    add = find_closest(reduced_vcf, cluster[1]) #find closest isolate
    cluster = c(cluster, add) #add next closest isolate to cluster
    if(in_list(sort(cluster),set_clusters)){ #if cluster is already picked, break out of loop
      break
    }
    clusters[[i]] <- cluster #save cluster
    num_unique_vars[i] = num_unique_variants(vcf, cluster, basemat) #calculate number of unique variants, store, etc.
    #print(num_unique_vars)
    
  }
  #pick clusters that are local maximia for number of unique variants and is above the min num of unique vars
  above_min = which(num_unique_vars >= min) #clusters with the minimum number of unique snps or more
  #add 2nd value to beginning and 2nd to last value to end of num_uniuqe_vars so first and last clusters can be local maxima
  num_unique_vars=c(num_unique_vars[2],num_unique_vars,num_unique_vars[length(num_unique_vars)-1]) 
  #find local maxima
  difs = diff(sign(diff(num_unique_vars))) 
  difs = which(difs == -2 | difs == -1)
  r = intersect(difs, above_min) #find clusters which are local maxima and fall above the min unique snps cutoff
  return(clusters[r]) #return clusters that fit criteria
}


in_list = function(vector, list){
  if(length(list) == 0){return(F)}
  for(i in 1:length(list)){
    x = list[[i]]
    if(is.null(x)){
      return(F)
    }
    else if(length(x) != length(vector)){
      next
    }
    else if(sum(x==vector) == length(vector)){
      return(T)
    }
  }
}

num_unique_variants = function(vcf, cluster, basemat){
  # calculates how many unique variants are in a cluster
  # inputs: 
  # vcf: matrix with columns of genomes and rows of variant positions
  # cluster: character vector with genome names or numerical vector with index of genome in vcf
  
  return(length(unique_variants(vcf, cluster, basemat)))
}

unique_variants = function(vcf, cluster, basemat){
  # identifies unique variants in a cluster
  # inputs: 
  # vcf: matrix with columns of genomes and rows of variant positions
  # cluster: character vector with genome names or numerical vector with index of genome in vcf
  
  #clean vcf
  if("ALT" %in% colnames(vcf)){
    vcf = vcf[!(grepl("N", vcf$ALT)),]
    vcf = vcf[!(grepl("N", vcf$REF)),]
    rownames(vcf) = vcf$POS
    vcf = vcf[,10:ncol(vcf)]
  }
  
  right_n = apply(basemat, 1, function(x){
    length(cluster) %in% x
  })
  vcf = vcf[which(right_n == T),]
  if(is.vector(vcf)){
    vcf = t(as.matrix(vcf))
  }
  rownames(vcf) = which(right_n == T)
  # if cluster is a character vector, convert to numerical
  if(is.character(cluster)){
    cluster = which(colnames(vcf) %in% cluster)
  }
  #calculate number of unique variants
  
  unique_vars = apply(vcf,1, function(x){ #for each variant in vcf
    clusterx = x[cluster] #bases in cluster
    non_cluster = x[setdiff(1:(length(x)),cluster)] #bases not in cluster
    if(length(unique(clusterx)) > 1){ #if more than one base in cluster
      return(FALSE) 
    }
    if(length(non_cluster)>0){ #if cluster is not the whole set in vcf
      return(length(intersect(clusterx,non_cluster)) == 0) #TRUE if bases in cluster don't intersect with bases in non-cluster
    }
    else{ #if all genomes are in cluster
      return(TRUE)  
    }
  })
  return(names(unique_vars)[which(unique_vars==TRUE)])
}


base_matrix = function(vcf){
  #calculate number of A's T's G's and C's for each variant
  bases = c("A", "T", "G", "C")
  basemat = matrix(nrow = nrow(vcf), ncol = 4)
  colnames(basemat) = bases
  basemat = sapply(bases, function(b){
    temp = ifelse(vcf == b, 1, 0)
    return(rowSums(temp))
  })
  return(basemat)
}

find_closest = function(vcf, genome){
  # finds closest genome to a given genome in vcf file
  # inputs: 
  # vcf: matrix with columns of genomes and rows of variant positions
  # cluster: character genome name or numerical with index of genome in vcf
  
  # if genome is a character vector, convert to numerical
  if(is.character(genome)){
    genome = which(colnames(vcf)==genome)
  }
  
  if(ncol(vcf) == 2){ #if only one other genome in vcf, return other genome name
    return(colnames(vcf)[-genome])
  }
  
  distance = apply(vcf[,-genome], 2, function(x){ #get distance of each genome to given genome
    return(sum(x != vcf[,genome])) #return number of bases that are different from the given genome
  })
  
  return(names(which(distance == min(distance)))[1]) #return one genome that has the minimum distance to the given genome
}

two_closest = function(vcf){
  # finds the two most related genomes in a vcf file
  # inputs: 
  # vcf: matrix with columns of genomes and rows of variant positions
  
  out_mat = matrix(ncol = ncol(vcf), nrow = ncol(vcf)) #initialize distance matrix- square with nrow, ncol = number of genomes
  c_list = 2:(nrow(out_mat)) #initialize columns to go through for each row
  for(r in 1:nrow(out_mat)){ #for each row in the distance matrix
    for(c in c_list){ #go through non-redundant columns
      out_mat[r,c] = sum(vcf[,r]==vcf[,c]) #set cell in distance matrix to be equal to the number of shared variants between the row and column genome
    }
    c_list = c_list[-1] #reduce which columns you go through for each row to avoid redundancy
  }
  max = which(out_mat == max(out_mat, na.rm = TRUE), arr.ind = TRUE)[1,] #find index of one cell with most shared variants
  return(colnames(vcf)[max]) #return genome names
}

reduce_vcf = function(vcf, cluster){
  # produce vcf matrix with only the variants that are shared by all of the genomes in a cluster
  # and with only the first genome in the cluster and all genomes not in cluster
  # inputs: 
  # vcf: matrix with columns of genomes and rows of variant positions
  # cluster: character genome name or numerical with index of genome in vcf
  cluster_vcf = vcf[,cluster] #subset of vcf with just genomes in cluster
  cluster_vcf_rows = apply(cluster_vcf, 1, function(x){ #for each variant
    return(length(unique(x)) == 1) #return if the same for all genomes in cluster
  })
  cluster_vcf_rows = which(cluster_vcf_rows == TRUE) 
  return_vcf = vcf[cluster_vcf_rows,] #vcf with all genomes and only variants that are all equal in cluster
  col_remove = cluster[-1] 
  col_remove = which(colnames(vcf) %in% col_remove)
  return_vcf = return_vcf[,-col_remove] #keep only one genome from cluster in the vcf
  return(return_vcf)
  
}

#####RECOMBINATION PERMUTATION TEST#######
recombination_events_per_gene = function(recombinations, genes, fullmsapos=NULL){
  #function to determine how many recombination events overlap with each gene
  #recombinations is data frame with 2 columns, one for start position and one for end position
  #genes is data frame with 2 columns, one for start position and one for end position
  #fullmsapos if recombination events in terms of core genome
  
  #convert recombination start and end positions to be in concordance with reference genome
  if(!is.null(fullmsapos)){
    recombinations[,1] = fullmsapos[recombinations[,1],]
    recombinations[,2] = fullmsapos[recombinations[,2],]
  }
  
  #define column names
  colnames(recombinations) = c("start", "end")
  colnames(genes) = c("start", "end")
  
  #calculate number of recombination events whose start and end positions are not both less than 
  #the gene start position or greater than the gene end position
  num = sapply(1:nrow(genes), function(c){      
    a = which(!(recombinations$start < genes$start[c] & recombinations$end < genes$start[c]))
    b = which(!(recombinations$start > genes$end[c] & recombinations$end > genes$end[c]))
    return(length(intersect(a,b)))
  })
  #return numerical vector of recombination events corresponding to each gene
  return(num)
}


place_rec_events = function(recombinations, core_genome_length=1754413){
  #function to randomly place recombination events throughout the genome based on length and number of inputs
  #recombinations is data frame with 2 columns, one for start position and one for end position
  #genome length is the length of the reference genome
  #fullmsapos if recombination events in terms of core genome
  
  colnames(recombinations) = c("start", "end")
  #determine length of each recombination event
  recombinations$length = recombinations$end-recombinations$start
  
  #choose random start positions for each recombination event
  recombinations$new_start = sample(1:core_genome_length, nrow(recombinations), replace = T)
  #make new random recombination event same length as original
  recombinations$new_end = recombinations$new_start + recombinations$length
  #re-calculate length
  recombinations$new_length = recombinations$new_end - recombinations$new_start
  
  #re-do new recombination events whose end is beyond the end of the reference genome
  no_end = which(recombinations$new_end > core_genome_length)
  for(i in no_end){
    rec_length = recombinations$length[i]
    start = sample(1:(core_genome_length-rec_length), 1)
    recombinations$new_start[i] = start
    recombinations$new_end[i] = start+rec_length
    recombinations$new_length = recombinations$new_end - recombinations$new_start
  }
  
  #return new recombination start and end positions
  return(recombinations[c("new_start", "new_end")])
  
}


rec_events_per_gene_permute = function(n, recombinations, genes, core_genome_length=1754413, fullmsapos=NULL){
  #function randomly places recombination events and calculates how many recombination events 
  #overlap with each gene n number of times
  #output is a matrix where each column corresponds to the inputted genes, and the rows are
  #the number of recombination events that overlap with that gene in each permutation
  
  #n = number of permutations
  #recombinations is data frame with 2 columns, one for start position and one for end position
  #genes is data frame with 2 columns, one for start position and one for end position
  #genome length is the length of the reference genome
  #fullmsapos if recombination events in terms of core genome
  
  #initialize output matrix
  output = matrix(nrow = n, ncol = nrow(genes))
  #permute n times
  for(i in 1:n){
    #print(i) #counter to keep track of which permutation is running
    write(i, file = "/nfs/esnitkin/Project_Cdiff/Analysis/targeted_sequencing_paper/2018-06-15_RM_variants_with_Ns_and_dashes/clades_1_and_2/counter", append = F)
    #randomly place recombination events
    random_rec_events = place_rec_events(recombinations, core_genome_length)
    #calculate number of overlapping recombination events for each gene
    random_rec_events_per_gene = recombination_events_per_gene(random_rec_events, genes, fullmsapos)
    #write to output matrix
    output[i,] = random_rec_events_per_gene
  }
  return(output)
}




rec_events_per_gene_pval = function(n, recombinations, genes, core_genome_length=1754413, fullmsapos=NULL){
  #function calculates empirical p-values for permutation test described in rec_events_per_gene_permute
  
  #n = number of permutations
  #recombinations is data frame with 2 columns, one for start position and one for end position
  #genes is data frame with 2 columns, one for start position and one for end position
  #genome length is the length of the reference genome
  #fullmsapos if recombination events in terms of core genome
  
  #calculate observed number of recombination events per gene
  observed_rec_events = recombination_events_per_gene(recombinations, genes, fullmsapos)
  #permute recombination events in random places
  permuted_rec_events = rec_events_per_gene_permute(n, recombinations, genes, core_genome_length, fullmsapos)
  #calculate p-value for each gene that gives probability of observing the number of observed recombination 
  #events if recombination events occured randomly throughout the genome
  pval = sapply(1:length(observed_rec_events), function(i){
    return((1+sum(permuted_rec_events[,i]>=observed_rec_events[i]))/(n+1))
  })
  out = rbind(permuted_rec_events, pval)
  return(out)
}

cl_list_to_clustervar = function(cl_list){
  clustervar = sapply(1:length(cl_list), function(x){
    return(rep(x,length(cl_list[[x]])))
  })
  clustervar = unlist(clustervar)
  names(clustervar) = unlist(cl_list)
  return(clustervar)
}



#######RECOMBINATION EVENTS PER BASE#################
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}


read.gff <- function(gubbins_gff_file){
  gff = gffRead(gubbins_gff_file)
  gff$node <- getAttributeField(gff$attributes, "node")
  gff$neg_log_likelihood <- getAttributeField(gff$attributes, "neg_log_likelihood")
  gff$taxa <- getAttributeField(gff$attributes, "taxa")
  gff$snp_count <- getAttributeField(gff$attributes, "snp_count")
  gff$snp_count <- gsub('"', "", gff$snp_count )
  gff$snp_count <- as.numeric(gff$snp_count)
  gff$length <- gff$end - gff$start
  return(gff)
  
}

rec_events_per_base <- function(gff_file, fullmsapos_file, genome_length=4290252){
  #Calculates how many recombination events overlap with each base in the reference genome
  #gff_file is the *.recombination_predictions.gff file from gubbins output

  #Read in gff and fullmsapos files
  gff = read.gff(gff_file)
  fullmsapos = read.table(fullmsapos_file)
  
  gff$ref_start = fullmsapos[gff$start,]
  gff$ref_end = fullmsapos[gff$end,]
  
  #Variable with all of the positions in each recombination event
  pos_rec_events = apply(gff, 1, function(x){
    return(x["ref_start"]: x["ref_end"])
  })
  
  #Combines positions of rec events into one long vector with each position in each event included
  #If a position is in 2 rec events, it will be in the list twice
  pos_rec_events = as.numeric(unlist(pos_rec_events))
  
  #Calculates how many times each position is in the vector
  pos_rec_events_table = table(pos_rec_events)
  
  #Add postitions without any recombination events to table
  w = which(!(0:genome_length %in% names(pos_rec_events_table)))
  w_table = table(w)
  w_table[0:length(w)] = 0
  pos_rec_events_table_with_zeros = c(w_table, pos_rec_events_table)
  pos_rec_events_table_with_zeros = pos_rec_events_table_with_zeros[order(as.numeric(names(pos_rec_events_table_with_zeros)))]
  
  return(pos_rec_events_table_with_zeros)
}
