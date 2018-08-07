#FUNCTIONS TO PROVIDE SUMMARIES OF GENOMIC DATA SETS
# this was called genomic_summary_v2_4.R

library(phytools)
library(ape)
library(geiger)
#v2_4
#minor fix in "collapse_pure_subtrees" function
#v2_3
# - fixed mixed_facility_subtree_4 and facility_clustering method 4
# - Added functions: 
#   1) collapse_pure_subtrees (collapses pure subtrees based on facility) 
#   2) get_pure_subtrees_v2 (returns all pure subtrees)
#v2
#	- Added new randomizations for calculating inter-facility clustering
#v2.1 
#2017-08-17: MODIFIED "pure_subtree_rand" BY MERGING PURE SUBTREES INTO SINGLE TIPS, AND RANDOMIZE ALL TIPS.
# - in case there are no pure subtrees, randomize all tips directly.
#2017-08-31: MODIFIED "pure_subtree_rand" BY PUTTING "subtree_tips" (IF ONLY ONE SUBTREE) INTO A LIST

#function intra_vs_inter_patient_variation
#
#Takes as input: 1) snp_dist    - A matrix of inter-sample SNP distances
#		 2) samples	- A 2 column matrix, where the first column is a sample ID and the second is patient ID
#		 3) prefix	- 
#
#Produces histogram of the number of variants within and between patients
intra_vs_inter_patient_variation = function(snp_dist, samples, prefix)
{
  
  #GET COMMON SAMPLE
  common_ids = intersect(samples[,1], row.names(snp_dist))
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[samples[,1] %in% common_ids, 1];
  sample_pts = samples[samples[,1] %in% common_ids, 2];
  
  snp_dist = snp_dist[sample_ids, sample_ids]
  
  #CALCULATE PAIRWISE SAMPLE DISTANCES
  snp_dist_1kSNPmax = snp_dist;
  snp_dist_1kSNPmax[snp_dist_1kSNPmax > 1000] = 1000;
  
  
  #GET PAIRS OF ISOLATES FROM THE SAME PATIENT
  unique_pts = unique(sample_pts);
  
  intra_pt_pairs = sapply(unique_pts, FUN = function(x){
    if(length(which(sample_pts == x)) > 1)
    {
      
      pt_samples = which(sample_pts == x);
      sample_pairs = combn(pt_samples, 2);
      
      sample_dist = snp_dist_1kSNPmax[t(sample_pairs)]
      
    }#end if
    
  });
  
  
  intra_pt_dist = unlist(intra_pt_pairs);
  
  print(paste(sum(unlist(lapply(intra_pt_pairs, length)) > 0), " patients have more than 1 sample sequenced." ,sep = ""))
  
  
  #PLOT HISTOGRAMS
  if(sum(unlist(lapply(intra_pt_pairs, length)) > 0) > 0)
  {
    dir.create("figures_no_postclust/",showWarnings = F)
    file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_intra_vs_inter_patient_pairwise_SNP_dist_histogram.pdf', sep = "");
    pdf(file);
    
    par(mfrow = c(2,1))
    hist(snp_dist_1kSNPmax[lower.tri(snp_dist_1kSNPmax)], col="red" , breaks = 100, main = "SNP distance among all pairs of isolates", xlab = "Distance (up to 1000)", xlim = c(0,max(snp_dist_1kSNPmax)))
    hist(intra_pt_dist, col="blue", breaks = 100, main = "SNP distance among isolates from the same patient", xlab = "Distance (up to 1000)", xlim = c(0,max(snp_dist_1kSNPmax)))
    
    dev.off();
    
  }#end if
  
}#end intra_vs_inter_patient_variation


#function intra_vs_inter_facility_variation
#
#Takes as input: 1) snp_dist	- A matrix of inter-sample SNP distances
#		 2) samples	- A 2 column matrix, where the first column is a sample ID and the second is the facility
#		 3) prefix	- A prefix for output files
#
#Produces histogram of the number of variants within and between patients
intra_vs_inter_facility_variation = function(snp_dist, samples, prefix)
{
  
  #GET COMMON SAMPLE
  common_ids = intersect(samples[,1], row.names(snp_dist))
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[samples[,1] %in% common_ids, 1];
  sample_facil = samples[samples[,1] %in% common_ids, 2];
  
  #CALCULATE PAIRWISE SAMPLE DISTANCES
  snp_dist = snp_dist[sample_ids, sample_ids]
  
  snp_dist_1kSNPmax = snp_dist;
  snp_dist_1kSNPmax[snp_dist_1kSNPmax > 1000] = 1000;
  
  
  #GET PAIRS OF ISOLATES FROM THE SAME PATIENT
  unique_facil = unique(sample_facil);
  
  intra_facil_pairs = sapply(unique_facil, FUN = function(x){
    if(length(which(sample_facil == x)) > 1)
    {
      
      facil_samples = which(sample_facil == x);
      sample_pairs = combn(facil_samples, 2);
      
      sample_dist = snp_dist_1kSNPmax[t(sample_pairs)]
      
    }#end if
    
  });
  
  
  intra_facil_dist = unlist(intra_facil_pairs);
  
  
  #PLOT HISTOGRAMS
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_intra_vs_inter_facility_pairwise_SNP_dist_histogram.pdf', sep = "");
  pdf(file);
  
  par(mfrow = c(2,1))
  hist(snp_dist_1kSNPmax[lower.tri(snp_dist_1kSNPmax)], col="red" , breaks = 100, main = "SNP distance among all pairs of isolates", xlab = "Distance (up to 1000)")
  hist(intra_facil_dist, col="blue", breaks = 100, main = "SNP distance among isolates from the same facility", xlab = "Distance (up to 1000)")
  
  dev.off();
  
}#end intra_vs_inter_facility_variation



#function intra_facility_variation_hist
#
#Takes as input: 1) snp_dist	- A matrix of inter-sample SNP distances
#		 2) samples	- A 2 column matrix, where the first column is a sample ID and the second is the facility
#		 3) prefix	- A prefix for output files
#		 4) plot_dim	- A two element vector indicating the dimensions of the multi-panel plot
#Produces histogram of the number of variants within and between patients
intra_facility_variation_hist = function(snp_dist, samples, prefix, plot_dim)
{
  
  #GET COMMON SAMPLE
  common_ids = intersect(samples[,1], row.names(snp_dist))
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[samples[,1] %in% common_ids, 1];
  sample_facil = samples[samples[,1] %in% common_ids, 2];
  
  #CALCULATE PAIRWISE SAMPLE DISTANCES
  snp_dist = snp_dist[sample_ids, sample_ids]
  
  snp_dist_1kSNPmax = snp_dist;
  snp_dist_1kSNPmax[snp_dist_1kSNPmax > 1000] = 1000;
  
  
  #GET PAIRS OF ISOLATES FROM THE SAME PATIENT
  unique_facil = unique(sample_facil);
  
  intra_facil_pairs = sapply(unique_facil, FUN = function(x){
    if(length(which(sample_facil == x)) > 1)
    {
      
      facil_samples = which(sample_facil == x);
      sample_pairs = combn(facil_samples, 2);
      
      sample_dist = snp_dist_1kSNPmax[t(sample_pairs)]
      
    }#end if
    
  });
  
  
  intra_facil_dist = unlist(intra_facil_pairs);
  
  
  #PLOT HISTOGRAMS
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_intra_facility_pairwise_SNP_dist_histogram.pdf', sep = "");
  pdf(file);
  
  par(mfrow = plot_dim)
  
  for(f in unique_facil)
  {
    
    hist(intra_facil_pairs[[f]], col="blue", breaks = 100, xlab = "Distance (up to 1000)", main = f, xlim = c(0, max(snp_dist_1kSNPmax)))
    
  }#end for
  
  dev.off();
  
}#end intra_facility_variation_hist


#function inter_facility_variation_hist
#
#Takes as input: 1) snp_dist	- A matrix of inter-sample SNP distances
#		 2) samples	- A 2 column matrix, where the first column is a sample ID and the second is the facility
#		 3) prefix	- A prefix for output files
#Produces histogram of pariwise genentic distance among isolates from all pairs of facilities
inter_facility_variation_hist = function(snp_dist, samples, prefix, plot_dim)
{
  
  #GET COMMON SAMPLE
  common_ids = intersect(samples[,1], row.names(snp_dist))
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[samples[,1] %in% common_ids, 1];
  sample_facil = samples[samples[,1] %in% common_ids, 2];
  
  #CALCULATE PAIRWISE SAMPLE DISTANCES
  snp_dist = snp_dist[sample_ids, sample_ids]
  
  snp_dist_1kSNPmax = snp_dist;
  snp_dist_1kSNPmax[snp_dist_1kSNPmax > 1000] = 1000;
  
  
  #GET PAIRS OF ISOLATES FROM THE SAME PATIENT
  unique_facil = unique(sample_facil);
  facil_pairs = t(combn(unique_facil,2))
  facil_pairs = rbind(facil_pairs, cbind(unique_facil, unique_facil))
  
  inter_facil_pairs = sapply(1:nrow(facil_pairs), FUN = function(x){
    if(length(which(sample_facil == facil_pairs[x,1])) > 1 && length(which(sample_facil == facil_pairs[x,2])) > 1)
    {
      
      sample_dist = snp_dist_1kSNPmax[sample_facil == facil_pairs[x,1], sample_facil == facil_pairs[x,2]]
      
    }#end if
    
  });
  
  
  
  #PLOT HISTOGRAMS
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_inter_facility_pairwise_SNP_dist_histogram.pdf', sep = "");
  pdf(file, height = 30, width = 30);
  
  par(mfrow = c(length(unique_facil), length(unique_facil)))
  
  for(p in 1:nrow(facil_pairs))
  {
    
    plot_x = which(unique_facil ==  facil_pairs[p,1])				
    plot_y = which(unique_facil ==  facil_pairs[p,2])				
    
    par(mfg = c(plot_x, plot_y))
    
    hist(inter_facil_pairs[[p]], col="blue", breaks = 100, xlab = "Distance (up to 1000)", main = paste(facil_pairs[p,1], ' vs ', facil_pairs[p,2], sep = ""), xlim = c(0, max(snp_dist_1kSNPmax)))
    
  }#end for
  
  dev.off();
  
}#end inter_facility_variation_hist



#function facility_fst
#
#Takes as input: 1) snp_dist	- A matrix of inter-sample SNP distances
#		 2) samples	- A 2 column matrix, where the first column is a patient ID and the second is the facility
#		 3) facils	- A list of facilities
#
#Computes a matrix of fst for each pair of facilities
facility_fst <- function(snp_dist, samples, facil)
{
  
  #GET COMMON SAMPLE
  common_ids = intersect(samples[,1], row.names(snp_dist))
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[samples[,1] %in% common_ids, 1];
  sample_facil = samples[samples[,1] %in% common_ids, 2];
  
  snp_dist = snp_dist[sample_ids, sample_ids]
  
  
  #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
  facil_dist = matrix(0, ncol = length(facil), nrow = length(facil), dimnames = list(facil, facil))
  
  for(f1 in unique(sample_facil))
  {
    
    for(f2 in unique(sample_facil))
    {
      
      facil_dist[f1,f2] = mean(snp_dist[sample_facil == f1, sample_facil == f2])
      
    }#end for
    
    
  }#end for
  
  return(facil_dist);
  
}#end facility_fst


#function facility_fst_cluster
#
#Takes as input: 1) snp_dist    - A matrix of inter-sample SNP distances
#                2) samples     - A 2 column matrix, where the first column is a patient ID and the second is the facility
#                3) facils      - A list of facilities
#                4) clusters    - A list of clusters to which samples belong to, named by sample ID
#
#Computes a matrix of fst for each pair of facilities, where intra/inter facility distances are restricted to samples within a cluster 
facility_fst_cluster <- function(snp_dist, samples, facil, clusters)
{
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[,1];
  sample_facil = samples[,2];
  
  snp_dist = snp_dist[sample_ids, sample_ids]
  clusters = clusters[sample_ids];
  
  clusters_gt10 = names(table(clusters))[table(clusters) >= 10]
  
  #CALCULATE INTRA- AND INTER-FACILITY DISTANCE
  facil_dist = matrix(0, ncol = length(facil), nrow = length(facil), dimnames = list(facil, facil))
  facil_dist_total = matrix(0, ncol = length(facil), nrow = length(facil), dimnames = list(facil, facil))
  facil_dist_count = matrix(0, ncol = length(facil), nrow = length(facil), dimnames = list(facil, facil))
  
  for(f1 in unique(sample_facil))
  {
    
    for(f2 in unique(sample_facil))
    {
      
      #for(c in unique(clusters))
      for(c in clusters_gt10)
      {
        
        if((sum((sample_facil == f1 & clusters == c) | (sample_facil == f2 & clusters == c)) >= 2) &&
           (sum(sample_facil == f1 & clusters == c) >=1) && (sum(sample_facil == f2 & clusters == c) >=1))
        {
          
          facil_dist_total[f1,f2] = facil_dist_total[f1,f2] + sum(snp_dist[sample_facil == f1 & clusters == c, sample_facil == f2 & clusters == c])
          facil_dist_count[f1,f2] = facil_dist_count[f1,f2] + sum(sample_facil == f1 & clusters == c)*sum(sample_facil == f2 & clusters == c);
          #facil_dist_total[f1,f2] = facil_dist_total[f1,f2] + mean(snp_dist[sample_facil == f1 & clusters == c, sample_facil == f2 & clusters == c])
          #facil_dist_count[f1,f2] = facil_dist_count[f1,f2] + 1;
          
        }#end if
        
      }#end for
      
    }#end for
    
  }#end for
  
  facil_dist = facil_dist_total / facil_dist_count;
  facil_dist[is.nan(facil_dist)] = 0;
  
  
  return(facil_dist);
  
}#end facility_fst_cluster

#facility_subtree - Takes as input: 1) subtrees          - A list produced by subtrees, that includes all isolates of interest
#                                   2) faciity	         - A named vactor of labels by which pure clusters are defined
#
#Returns the number of putative intra-facility transmissions based on clustering of isolates from same facility
facility_subtree <- function(subtrees, facility)
{
  
  largest_st = numeric();
  
  #GET LARGEST SUB-TREE FOR EACH ISOLATE
  for(i in as.character(names(facility)))
  {
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVE ONLY A SINGLE FACILITY LABEL
    largest_st[i] = max(sapply(subtrees, FUN = function(st){ if(sum(i %in% st$tip.label) > 0   && 
                                                                length(unique(facility[as.character(intersect(st$tip.label, names(facility)))])) == 1 )
    {
      length(intersect(names(facility), st$tip.label))
    }else{
      0
    }
    }));
    
  }#end for
  
  #COUNT NUMBER OF PUTATIVE TRANSMISSIONS (SUBTRACT 1 FROM EACH SUB-TREE TO ACCOUNT FOR INTRODUCTION INTO FACILITY)
  facil_tns = structure(rep(0, length(unique(facility))) , names = sort(unique(facility)));
  
  for(f in unique(facility))
  {
    
    st_count = table(largest_st[facility == f])
    print(st_count)
    
    for(c in names(st_count)[as.numeric(names(st_count)) > 0])
    {
      
      facil_tns[f] = facil_tns[f] + (st_count[c]/as.numeric(c))*(as.numeric(c)-1);
      
    }#end for
    
  }#end for
  
  return(facil_tns)
  
}#end facility_subtree


#patient_subtree - Takes as input: 1) subtrees    - A list produced by subtrees, that includes all isolates of interest
#                                  2) patient     - A named vactor of patient IDs
#
#Returns a set of patients, with a single representative of patient clusters selected
patient_subtree <- function(subtrees, patient)
{
  
  largest_st = structure(rep(0,length(patient)), names = as.character(names(patient)));
  largest_st_i = structure(rep(0,length(patient)), names = as.character(names(patient)));
  
  #GET LARGEST SUB-TREE FOR EACH ISOLATE
  for(i in as.character(names(patient)))
  {
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVE ONLY A SINGLE FACILITY LABEL
    largest_st[i] = max(sapply(subtrees, 
                               FUN = function(st){ 
                                 if(sum(i %in% st$tip.label) > 0   && 
                                    length(unique(patient[as.character(intersect(st$tip.label, names(patient)))])) == 1  &&
                                    length(intersect(names(patient), st$tip.label)) == length(st$tip.label))
                                   
                                 {
                                   length(intersect(names(patient), st$tip.label))
                                 }else{
                                   0
                                 }
                               }));
    
    largest_st_i[i] = which.max(sapply(subtrees,  
                                       FUN = function(st){ 
                                         if(sum(i %in% st$tip.label) > 0   && 
                                            length(unique(patient[as.character(intersect(st$tip.label, names(patient)))])) == 1 &&
                                            length(intersect(names(patient), st$tip.label)) == length(st$tip.label))
                                           
                                         {
                                           length(intersect(names(patient), st$tip.label))
                                         }else{
                                           0
                                         }
                                       }));
    
    
  }#end for
  
  
  #GET SET OF ISOLATES, WITH A SINGLE REPRESENTATIVE TAKEN FOR PATIENTS WHOSE ISOLATES CLUSTER
  pt_subset = names(largest_st)[largest_st == 0]
  
  for(st in setdiff(unique(largest_st_i), 1))
  {
    pt_subset = c(pt_subset, subtrees[[st]]$tip.label[1])
    
  }#end for
  
  
  print(length(largest_st))
  print(sum(largest_st == 0))
  print(length(setdiff(unique(largest_st_i), 1)))
  print(length(pt_subset))
  
  return(pt_subset)
  
}#end patient_subtree


#function intra_pt_var_hist
#
#Takes as input: 1) var_aln     - A var object
#                2) samples     - A 2 column matrix, where the first column is a patient ID and the second is sample ID
#
#Produces histograms of variants across the genomes of samples from the same patient to look for evidence of recombination/sequence errors
intra_pt_var_hist <- function(var_aln, samples)
{
  #GET COMMON SAMPLE
  common_ids = intersect(samples[,1], row.names(var_aln))
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[samples[,1] %in% common_ids, 1];
  sample_facil = samples[samples[,1] %in% common_ids, 2];
  
  var_aln = var_aln[sample_ids,]
  
  #CALCULATE PAIRWISE SAMPLE DISTANCES
  var_dist = dist.dna(var_aln, model = "raw");
  
  snp_dist = as.matrix(var_dist) * ncol(var_aln);
  
  snp_dist_1kSNPmax = snp_dist;
  snp_dist_1kSNPmax[snp_dist_1kSNPmax > 1000] = 1000;
  
  
  #GET PAIRS OF ISOLATES FROM THE SAME PATIENT
  unique_pts = unique(sample_pts);
  
  intra_pt_pairs = sapply(unique_pts, FUN = function(x){
    if(length(which(sample_pts == x)) > 1)
    {
      
      pt_samples = which(sample_pts == x);
      sample_pairs = combn(pt_samples, 2);
      
      sample_dist = snp_dist_1kSNPmax[t(sample_pairs)]
      
    }#end if
    
  });
  
  
  intra_pt_dist = unlist(intra_pt_pairs);
  
  print(paste(sum(unlist(lapply(intra_pt_pairs, length)) > 0), " patients have more than 1 sample sequenced." ,sep = ""))
  
  
  #PLOT INTRA-PATIENT HISTOGRAMS OF VARIANTS ACROSS THE GENOME
  intra_pt_dates = numeric();
  
  for(p in unique_pts)
  {
    
    #CHECK IF MULTIPLE SAMPLES FOR CURRENT PATIENT
    if(length(which(sample_pts == p)) > 1)
    {
      
      #GET SAMPLES ASSOCIATED WITH CURRENT PATIENT
      pt_samples = which(sample_pts == p);
      
      #PLOT HISTOGRAMS OF SNP DENSITY FOR EACH PAIR OF ISOLATES
      dir.create("figures_no_postclust/",showWarnings = F)
      file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_pt_', p, '_intra_patient_snp_dist.pdf', sep = "");
      pdf(file);
      par( mfrow = c(length(pt_samples), length(pt_samples)))
      par(oma = c (4, 4, 4, 0))
      
      sample_pairs = combn(pt_samples, 2);
      
      for(sp in 1:ncol(sample_pairs))
      {
        
        snps = which(dna[samples_i[sample_pairs[1,sp]] ,] != dna[samples_i[sample_pairs[2,sp]],]);
        
        
        par(mfg = c(which(pt_samples == sample_pairs[1,sp]), which(pt_samples == sample_pairs[2,sp])))
        par(mar = c(2,2,2,2));
        
        if(length(snps) > 0)
        {
          hist(snps, xlim = c(1, ncol(dna)), 10000, main = paste(sample_names_subset[sample_pairs[1,sp]], " vs ", sample_names_subset[sample_pairs[2,sp]]))
          
        }else{
          
          hist(c(0), xlim = c(1, ncol(dna)), 10000, main = paste(sample_names_subset[sample_pairs[1,sp]], " vs ", sample_names_subset[sample_pairs[2,sp]]))
          
        }#end if
        
        #PRINT OUT INFORMATION FOR INTRA-PATIENT SAMPLE PAIRS WITH HIGH SNP DISTANCE
        if(length(snps) > 50)
        {
          culture1 = evan_sample_key$ Culture.Result[evan_sample_key$Sample.ID == sample_names_subset[sample_pairs[1,sp]]];
          date1 =  evan_sample_key$Culture.Date[evan_sample_key$Sample.ID == sample_names_subset[sample_pairs[1,sp]]];
          
          culture2 = evan_sample_key$ Culture.Result[evan_sample_key$Sample.ID == sample_names_subset[sample_pairs[2,sp]]];
          date2 =  evan_sample_key$Culture.Date[evan_sample_key$Sample.ID == sample_names_subset[sample_pairs[2,sp]]];
          
          intra_pt_dates = c(intra_pt_dates, abs(as.numeric(as.character(date2)) - as.numeric(as.character(date1))))
          print(paste(p, sample_names_subset[sample_pairs[1,sp]],culture1,date1 ,sample_names_subset[sample_pairs[2,sp]], culture2, date2, length(snps) , sep = " , "))
          
          
        }#end if
        
      }#end for
      
      mtext(paste("Patient ", p), side = 3, outer = T, line = 2)
      dev.off()
      
      
    }#end if
    
  }#end for
  
}#end intra_pt_var_hist


#function snp_hist
#
#Takes as input: 1) dna         - A DNAaln object
#                2) snp_thresh  - A threshold to break clusters
#                3) prefix      - A prefix for output figures
#
#Shows distribution of variants across genome for all sequences and within each cluster
snp_hist <- function(dna, snp_thresh, prefix)
{
  
  #GET VARIABLE POSITIONS ACROSS ALL SEQUENCES AND PLOT DISTRIBUTION
  var_pos = apply(dna, 2, FUN = function(x)
  {
    
    sum(x != x[1] | x == 'N') > 0;
    
  })
  
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_variant_hist.pdf', sep = "");
  
  pdf(file)
  
  hist(which(var_pos), 5000, xlab = "Genome position", ylab = "Number of variants")
  
  dev.off();
  
  
  #GET SEQUENCE CLUSTERS
  snp_dist = matrix(outer(1:nrow(dna), 1:nrow(dna), FUN = Vectorize(function(x,y){sum(dna[x, var_pos] != dna[y, var_pos])})),
                    nrow = nrow(dna), ncol = nrow(dna), dimnames = list(row.names(dna), row.names(dna)));
  
  snp_clust = hclust(as.dist(snp_dist))
  
  clusters = cutree(snp_clust,h = snp_thresh)
  
  
  #PLOT VARIABLE POSITIONS WIHTIN EACH CLUSTER
  for(c in unique(clusters))
  {
    
    if(sum(clusters == c) > 1)
    {
      
      c_names = names(clusters)[clusters == c]
      
      c_var_pos = apply(dna[c_names,], 2, FUN = function(x)
      {
        
        #sum(x != x[1] || x == 'N') > 0; #CHANGE FROM || TO |
        sum(x != x[1] | x == 'N') > 0
        
      })
      
      dir.create("figures_no_postclust/",showWarnings = F)
      file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, 'cluster_', c, '_numMemb', length(c_names), '_variant_hist.pdf', sep = "");
      
      pdf(file)
      
      hist(which(c_var_pos), 5000, xlab = "Genome position", ylab = "Number of variants")
      
      dev.off();
      
    }#end if
    
  }#end
  
  
}#end snp_hist

#LOOK AT RECOMBINATION-FILTERED SNPS

#function pw_dist_filt
#
#Takes as input: 1) var_aln		- A DNA object representing a multiple sequence alignment
#
#Gets pairwise genetic distances among isolates
pw_dist = function(var_aln)
{
  
  #GET SNP DISTANCE MATRIX
  var_dist = dist.dna(var_aln[,filt_pos], model = "raw");
  
  snp_dist_mat = as.matrix(var_dist) * length(filt_pos);
  snp_dist_mat[snp_dist_mat == -Inf] = 0;
  
  return(snp_dist_mat)
  
}#end pw_dist


#function pw_dist_filt
#
#Takes as input: 1) var_aln		- A DNA object representing a multiple sequence alignment
#		 2) filt_pos		- A vector of positions to consider
#
#Gets pairwise genetic distances among isolates with positions filtered out
pw_dist_filt = function(var_aln, filt_pos)
{
  
  #GET SNP DISTANCE MATRIX
  var_dist = dist.dna(var_aln[,filt_pos], model = "raw");
  
  snp_dist_mat = as.matrix(var_dist) * length(filt_pos);
  snp_dist_mat[snp_dist_mat == -Inf] = 0;
  
  return(snp_dist_mat)
  
}#end pw_dist_filt

#function snp_hist_filt
#
#Takes as input: 1) dna		- A DNAaln object
#		 2) snp_thresh	- A threshold to break clusters
#		 3) prefix	- A prefix for output figures
#		 4) filt_snps	- A list of positions to filter out
#
#Shows distribution of variants across genome for all sequences and within each cluster
snp_hist_filt <- function(dna, snp_thresh, prefix, filt_snps)
{
  
  #GET VARIABLE POSITIONS ACROSS ALL SEQUENCES AND PLOT DISTRIBUTION
  var_pos = apply(dna, 2, FUN = function(x)
  {
    
    sum(x != x[1] | x == 'N') > 0;
    
  })	
  
  print(c(sum(var_pos), length(filt_snps)))
  
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_variant_hist.pdf', sep = "");
  
  pdf(file)
  
  hist(intersect(which(var_pos), filt_snps), 5000, xlab = "Genome position", ylab = "Number of variants")
  
  dev.off();
  
  
  #GET SEQUENCE CLUSTERS
  snp_dist = matrix(outer(1:nrow(dna), 1:nrow(dna), FUN = Vectorize(function(x,y){sum(dna[x, var_pos] != dna[y, var_pos])})),
                    nrow = nrow(dna), ncol = nrow(dna), dimnames = list(row.names(dna), row.names(dna)));
  
  snp_clust = hclust(as.dist(snp_dist))
  
  clusters = cutree(snp_clust,h = snp_thresh)
  
  
  #PLOT VARIABLE POSITIONS WIHTIN EACH CLUSTER
  for(c in unique(clusters))
  {
    
    if(sum(clusters == c) > 1)
    {
      
      c_names = names(clusters)[clusters == c]
      
      c_var_pos = apply(dna[c_names,], 2, FUN = function(x)
      {
        
        sum(x != x[1] | x == 'N') > 0;
        
      })
      
      dir.create("figures_no_postclust/",showWarnings = F)
      file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, 'cluster_', c, '_numMemb', length(c_names), '_variant_hist.pdf', sep = "");
      
      pdf(file)
      
      hist(intersect(which(c_var_pos), filt_snps), 5000, xlab = "Genome position", ylab = "Number of variants")
      
      dev.off();
      
    }#end if
    
  }#end
  
  
}#end snp_hist_filt


#function tree_compare
#
#Takes as input: 1) tree1       - A phylogenetic tree
#                2) tree2       - A phylogenetic tree
#                3) prefix      - A prefix for naming outout plots
#
#Produces summary plots comparing the trees and distance matrices for the two trees
tree_compare <- function(tree1, tree2, prefix)
{
  
  #GET DISTANCE MATRICES
  d_mat1 = cophenetic.phylo(tree1);
  d_mat2 = cophenetic.phylo(tree2);
  
  common_seqs = intersect(row.names(d_mat1), row.names(d_mat2));
  inds = t(combn(as.character(common_seqs), 2))
  
  
  #PLOT COMPARISON OF MATRICES
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_distance_mat_scatter.pdf', sep = "");
  
  pdf(file)
  
  cor = cor.test(d_mat1[inds], d_mat2[inds])
  
  plot(d_mat1[inds], d_mat2[inds], xlab = "Distance mat 1", ylab = "Distance matrix 2")
  text(max(d_mat1[inds])*0.2, max(d_mat2[inds])*0.8, paste('R = ', round(cor$estimate,2), '(p =', round(cor$p.value,3), ')', sep = ""));
  
  dev.off();
  
  
}#end tree_compare


#function plot_cluster_trees
#
#Takes as input: 1) tree		- A phylogeny including all isolates
#		 2) clusters		- A vector the length tree tips indicating a cluster number
#		 3) tip_legend		- A vector of categorical assignments named by tree labels
#		 4) tip_colors		- A named vector of colors corresponding to categories in legend
#		 5) prefix		- A prefix for output files
#		 6) tree_labels		- An optional vector named by tip labels, with display tip labels
#
#Plots trees for isolate sub-clusters and returns sub-clusters
plot_cluster_trees = function(tree, clusters, tip_legend, tip_colors, prefix, tree_labels = NULL)
{
  #GO THROUGH EACH CLUSTER AND PLOT SUBTREE
  for(c in unique(clusters))
  {
    if(sum(clusters == c) > 1)
    {
      #GET CURRENT CLUSTER
      #cluster_mrca = getMRCA(reroot(tree , which.max(rowMeans(cophenetic(tree)))), names(clusters)[clusters == c]);
      #cluster_mrca = getMRCA(tree, names(clusters)[clusters == c])
      #sub_clade = extract.clade(tree, cluster_mrca); 
      sub_clade = drop.tip(tree, setdiff(tree$tip.label, names(clusters)[clusters == c]))
      
      tip_legend_subset = tip_legend[sub_clade$tip.label]
      
      print(paste(c, ": ", length(sub_clade$tip.label), sep = ""));	
      
      #RENAME TIPS
      if(!is.null(tree_labels))
      {
        
        sub_clade$tip.label = tree_labels[sub_clade$tip.label ];
        
      }#end if
      
      
      #PLOT TREE
      dir.create("figures_no_postclust/",showWarnings = F)
      file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_cluster', c, '_convertTree.pdf', sep = "");
      pdf(file, width = 6, height = 6)
      
      scaling_ratio =(length(tree$tip.label)-length(sub_clade$tip.label))/length(tree$tip.label);
      
      plot(sub_clade, type = "phylogram", label.offset = scaling_ratio*0.000005, no.margin = TRUE, cex = scaling_ratio)
      
      tiplabels(pie = to.matrix(tip_legend_subset,names(tip_colors)),piecol = tip_colors, cex = scaling_ratio*0.5)
      
      legend('topright', legend = names(tip_colors), col = "black", pt.bg = tip_colors, pch = 21)
      
      dev.off();
      
    }#end if
    
  }#end for
  
  
  #PLOT FULL TREE
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_full_convertTree.pdf', sep = "");
  pdf(file, width = 6, height = 6)
  
  tip_legend_subset = tip_legend[tree$tip.label];
  
  plot(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE)
  tiplabels(pie = to.matrix(tip_legend_subset,names(tip_colors)),piecol = tip_colors, cex = 0.2)
  legend('topright', legend = names(tip_colors), col = "black", pt.bg = tip_colors, pch = 21)
  
  dev.off();
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_full_convertTree_noEL.pdf', sep = "");
  pdf(file, width = 6, height = 6)
  
  tree$edge.length = log(tree$edge.length + 1)
  plot(tree, type = "fan", show.tip.label = FALSE, no.margin = TRUE)
  tiplabels(pie = to.matrix(tip_legend_subset,names(tip_colors)),piecol = tip_colors, cex = 0.2)
  legend('topright', legend = names(tip_colors), col = "black", pt.bg = tip_colors, pch = 21)
  
  dev.off();
  
}#end plot_cluster_trees


#function pw_dist_cluster
#
#Takes as input: 1) var_aln		- A DNA object representing a multiple sequence alignment
#		 2) isolate_legend	- A vector of categorical assignments named by tree labels
#		 3) isolate_colors	- A named vector of colors corresponding to categories in legend
#
#Plots heatmap of pairwise genetic distances among isolates
pw_dist_cluster = function(var_aln, isolate_legend, isolate_colors, prefix)
{
  
  #GET SNP DISTANCE MATRIX
  var_dist = dist.dna(var_aln, model = "raw");
  
  snp_dist_mat = as.matrix(var_dist) * ncol(var_aln);
  snp_dist_mat[snp_dist_mat == -Inf] = 0;
  
  
  #CREATE A LOG DISTANCE MATRIX FOR VISIBILITY
  log_snp_dist_mat = log2(snp_dist_mat);
  log_snp_dist_mat[log_snp_dist_mat == -Inf] = 0;
  
  #CREATE ROW COLUMN VECTOR
  rowCol = structure(isolate_colors[isolate_legend[names(row.names(log_snp_dist_mat))]], names = row.names(log_snp_dist_mat))	
  
  #PLOT HEATMAP
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), "_", prefix, '_log2_pw_genetic_distance_heatmap.pdf', sep = "");
  pdf(file, width = 6, height = 6)
  
  heatmap(log_snp_dist_mat, scale = "none", ColSideColors = rowCol, RowSideColors = rowCol, cexCol = min(1, 40/nrow(log_snp_dist_mat)), cexRow = min(1, 40/nrow(log_snp_dist_mat)))
  
  dev.off();
  
  return(snp_dist_mat)
  
}#end pw_dist_cluster


#function intra_pt_distance_scatter
#
#Takes as input: 1) var_aln	- A var object
#		 2) samples	- A 2 column matrix, where the first column is a patient ID and the second is sample ID
#
#Produces a scatter plot for each patient showing the pairwise distance between their isolates
intra_pt_distance_scatter = function(var_aln, samples)
{
  
  #GET COMMON SAMPLE
  common_ids = intersect(samples[,1], row.names(var_aln))
  
  #ASSIGN ARGUMENTS TO VARIABLES
  sample_ids = samples[samples[,1] %in% common_ids, 1];
  sample_pts = samples[samples[,1] %in% common_ids, 2];
  
  var_aln = var_aln[common_ids,]
  
  #CALCULATE PAIRWISE SAMPLE DISTANCES
  var_dist = dist.dna(var_aln, model = "raw");
  
  snp_dist = log2(as.matrix(var_dist) * ncol(var_aln));
  snp_dist[snp_dist == -Inf] = 0;
  
  
  #GET PAIRS OF ISOLATES FROM THE SAME PATIENT
  unique_pts = unique(sample_pts);
  
  intra_pt_isolates = sapply(unique_pts, FUN = function(x){
    if(length(which(sample_pts == x)) > 1)
    {
      
      pt_samples = which(sample_pts == x);
      
    }#end if
    
  });
  
  
  #PLOT DISTANCE SCATTER FOR EACH ISOLATE FROM A PATIENT WITH MULTIPLE ISOLATES
  dir.create("figures_no_postclust/",showWarnings = F)
  file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_log2_intra-patient_pw_SNV_distance_scatter.pdf', sep = "");
  pdf(file, width = 22, height = 7)
  
  count = 1;
  
  num_pt_isolates = unlist(lapply(intra_pt_isolates, length))
  for(pt in which(num_pt_isolates > 0))
  {
    
    #GET CURRENT PATIENT ISOLATES
    pt_inds = intra_pt_isolates[[pt]];
    
    
    #PLOT DISTANCE SCATTER FOR EACH PATIENT ISOLATE, HIHGLIGHTING OTHER ISOLATES FROM SAME PATIENT
    for(i in pt_inds)
    {
      #GET INDICIES FOR WITHIN AND BETWEEN PATIENTS
      intra_pt_i = setdiff(pt_inds, i);
      inter_pt_i = setdiff(1:nrow(snp_dist), pt_inds);
      
      plot(jitter(rep(count, length(intra_pt_i)), 0.2), (snp_dist[i, intra_pt_i]), xlim = c(0, sum(num_pt_isolates)+1), ylim = c(0, max(snp_dist)), axes = FALSE, pch = 17, col = "red", xlab = "", ylab = "");
      par(new = TRUE)
      plot(jitter(rep(count, length(inter_pt_i)), 0.2), (snp_dist[i, inter_pt_i]), xlim = c(0, sum(num_pt_isolates)+1), ylim = c(0,max(snp_dist)), axes = FALSE, pch = 21, col = "black", xlab = "", ylab = "");
      par(new = TRUE)
      
      #INCREMENT PATIENT COUNTER
      count = count + 1;
      
    }#end for
    
    
  }#end for
  
  
  axis(1, at = 1:sum(num_pt_isolates), labels = row.names(snp_dist)[unlist(intra_pt_isolates)], las = 2)
  axis(2, at = 0:max(snp_dist), las = 2);
  mtext("log2 SNV Distance", side=2, line=3, , las=3)
  
  dev.off();
  
}#end intra_pt_distance_scatter


########################### FUNCTIONS FOR TO DETERMINE EMPIRIC P-VALUE OF CLUSTERING OF ISOLATES FROM FACILITY PAIRS ##########################

#facility_clustering - Takes as input: 1) clusters      - A named vector indicating the clusters to which each isolate belongs
#                                      2) tree          - A phylogeny
#                                      3) facility      - A named vactor of facility labeled
#                                      4) num_perm      - Number of permutations on the tree to do
#                                      5) prefix        - A prefix for plotting
#				       6) facil_subset  - A subset of facilities for plotting all vs. all
#				       7) method	 - The method for calculating inter-facility p-values
#							   1 = Metric - # isolates from both facilities in mixed subtrees
#							       Rand   - Fix isolates from  first facility and randomize rest
#							       *Assumes reciprocal p-values will be combined
#							   2 = Metric - # isolates from first facilty
#							       Rand - Fix isolates from second facility
#							       *Assumes reciprocal p-values will be combined
#							   3 = Metric - Number of mixed subtrees 
#							       Rand - Swap pure subtrees
#							   4 = Metric - Maximal number of mixed subtrees (deconstruct mixed sub-trees)
#							       Rand - Swap pure subtrees  
#Returns the p-value associated with isolates from 2 facilities clustering more than expectation on the phylogeny
facility_clustering <- function(clusters, tree, facility, num_perm, prefix, facil_subset = NULL, method)
{
  
  facility = facility[tree$tip.label];
  
  facils = sort(unique(facility));
  
  facil_p_all = list();
  
  if(is.null(facil_subset))
  {
    facil_subset = facils;
    
  }#end if
  
  
  #GET CLUSTERS WITH AT LEAST 10 INDIVIDUALS ON THE TREE
  clusters = clusters[tree$tip.label];
  
  clusters_gt20 = names(table(clusters))[table(clusters) >= 0]
  
  
  #FOR EACH CLUSTER CALCULATE THE ASSOCIATION BETWEEN ISOLATES FROM EACH PAIR OF FACILITIES
  for(c in clusters_gt20)
  {
    
    facil_p = matrix(1, nrow = length(facils), ncol = length(facils),
                     dimnames = list(facils, facils))
    
    #GET TREE WITH ISOLATES IN CURRENT CLUSTER
    cluster_tree = drop.tip(tree, setdiff(tree$tip.label, names(clusters)[clusters == c]))
    
    
    #GET SUBTREES
    st = subtrees(cluster_tree)
    
    
    #CALCULATE NUMBER OF ISOLATES IN FROM FACILITY PAIRS GROUPED ON TREE IN REAL AND PERMUTED DATA
    facil_ints = list();
    
    if(method == 1)
    {
      
      #METHOD - COUNT TOTAL NUMBER OF ISOLATES FROM BOTH FACILITIES IN MIXED SUBTREES
      #RAND	- FIX ISOLATES FROM FIRST FACILITY AND RANDOMIZE THE REST
      
      facil_ints[[1]] = mixed_facility_subtree_1(st, facility[cluster_tree$tip.label])
      
      for(r in 1:num_perm)
      {
        print(r);
        
        facil_ints[[r+1]] = mixed_facility_subtree_1_rand(st, facility[cluster_tree$tip.label])
        
      }#end for
      
    }else if(method == 2){
      
      #METHOD - COUNT TOTAL NUMBER OF ISOLATES FROM FIRST FACILITY IN MIXED SUBTREES
      #RAND	- FIX ISOLATES FROM SECOND FACILITY
      facil_ints[[1]] = mixed_facility_subtree_2(st, facility[cluster_tree$tip.label])
      
      for(r in 1:num_perm)
      {
        print(r);
        
        facil_ints[[r+1]] = mixed_facility_subtree_2_rand(st, facility[cluster_tree$tip.label])
        
      }#end for
      
    }else if(method == 3){
      
      #METHOD - COUNT TOTAL NUMBER OF MIXED SUBTREES
      #RAND	- RANDOMIZE LOCATION OF PURE SUBTREES
      facil_ints[[1]] = mixed_facility_subtree_3(cluster_tree, st, facility[cluster_tree$tip.label])
      
      for(r in 1:num_perm)
      {
        print(r);
        cluster_tree_rand = pure_subtree_rand(cluster_tree, facility)
        cluster_tree_rand_st = subtrees(cluster_tree_rand)
        
        facil_ints[[r+1]] = mixed_facility_subtree_3(cluster_tree_rand, cluster_tree_rand_st, facility[cluster_tree_rand$tip.label])
        
      }#end for
      
    }else if(method == 4){
      
      #METHOD - COUNT MAXIMAL NUMBER OF MIXED SUBTREES (CONSIDERS EMBEDDED SUBTREES)
      #RAND	- RANDOMIZE LOCATION OF PURE SUBTREES
      #collapse first
      cluster_tree_pure = collapse_pure_subtrees(cluster_tree,facility[cluster_tree$tip.label])
      cluster_tree_pure_sts = subtrees(cluster_tree_pure)
      facil_ints[[1]] = mixed_facility_subtree_4(cluster_tree_pure, cluster_tree_pure_sts, facility[cluster_tree_pure$tip.label])
      
      for(r in 1:num_perm)
      {
        print(r);
        
        cluster_tree_rand = pure_subtree_rand(cluster_tree_pure, facility[cluster_tree_pure$tip.label])
        cluster_tree_rand_st = subtrees(cluster_tree_rand)
        
        facil_ints[[r+1]] = mixed_facility_subtree_4(cluster_tree_rand, cluster_tree_rand_st, facility[cluster_tree_rand$tip.label])
        
      }#end for
      
    }else if(method == 5){
      
      #METHOD - COUNT TRANSMISSIONS BASED ON ML ANCESTRAL RECONSTRUCTION
      #RAND	- RANDOMIZE LOCATION OF FACILITIES
      facil_ints[[1]] = get_facility_genet(cluster_tree, facility[cluster_tree$tip.label])
      
      for(r in 1:num_perm)
      {
        print(r);
        
        rand_facility = facility;
        rand_facility[names(facility)] = sample(facility, length(facility))
        
        facil_ints[[r+1]] = get_facility_genet(cluster_tree, rand_facility)
        
      }#end for
      
    }else if(method == 6){
      
      #METHOD - COUNT TRANSMISSIONS BASED ON ML ANCESTRAL RECONSTRUCTION
      #RAND	- RANDOMIZE LOCATION OF PURE SUBTREES
      facil_ints[[1]] = get_facility_genet(cluster_tree, facility[cluster_tree$tip.label])
      
      for(r in 1:num_perm)
      {
        print(r);
        
        cluster_tree_rand = pure_subtree_rand(cluster_tree, facility)
        cluster_tree_rand_st = subtrees(cluster_tree_rand)
        
        facil_ints[[r+1]] = get_facility_genet(cluster_tree_rand, facility[cluster_tree_rand$tip.label],facil_thresh = 0, BS_thresh = 0)
        
      }#end for
      
    }else{
      
      print(paste("Unkown method " , method, " selected!!!!", sep = ""))
      
      return(-1);
      
    }#end if
    
    
    #CALCULATE THE FRACTION OF PERMUTATIONS WITH MORE THAN THE OBSERVED NUMBER OF GROUPED ISOLATES FOR EACH FACILITY PAIR
    dir.create("figures_no_postclust/",showWarnings = F)
    file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_cluster_', c, '_method_', method, '_phylo_cluster_rand_hist.pdf', sep = "");
    pdf(file, height = 50, width = 50);
    
    par(mfrow = c(length(facil_subset), length(facil_subset)))
    
    for(i in 1:length(facil_subset))
    {
      
      for(j in 1:length(facil_subset))
      {
        
        counts = unlist(lapply(facil_ints, FUN = function(x){return(x[facil_subset[i], facil_subset[j]])}))
        
        facil_p[facil_subset[i], facil_subset[j]] = max(1/num_perm, sum(counts[1] < counts[2:(num_perm+1)]) / num_perm);
        
        h = hist(counts[2:(num_perm+1)], max(counts) + 1, plot = FALSE);
        
        plot(h, xlim = c(0,max(counts)+1), ylim = c(0, max(h$counts)+1), xlab = 'Tips',
             main = paste(facil_subset[i],' and ', facil_subset[j] ,sep = ""))
        
        par(new = TRUE);
        
        plot(counts[1], 0, cex = 2, col = 'red', pch = 18, axes = FALSE, xlab = "", ylab = "",
             xlim = c(0,max(counts)+1), ylim = c(0, max(h$counts)+1))
        
      }#end for
      
    }#end for
    
    dev.off();
    
    
    #PLOT RECIPROCAL INTERACTIONS AGAINST ONE ANOTHER
    dir.create("figures_no_postclust/",showWarnings = F)
    file = paste("figures_no_postclust/", format(Sys.time(), "%Y-%m-%d"), '_', prefix, '_cluster_', c, '_method_', method, '_recip_log10p_scatter.pdf', sep = "");
    pdf(file);
    
    plot(-1*log10(facil_p[t(combn(facil_subset,2))]) , -1*log10(facil_p[t(combn(facil_subset,2))[,2:1]]))
    
    dev.off()
    
    #SAVE VALUES TO RETURN
    facil_p_all[[c]] = facil_p;
    
  }#end for
  
  
  return(facil_p_all)
  
}#end facility_clustering


#FACILITY CLUSTERING METHOD 1
#mixed_facility_subtree_1 - Takes as input: 1) subtrees          - A list produced by subtrees, that includes all isolates of interest
#                                           2) facility          - A named vactor of facility labeled
#
#For each pair of facilities returns the number of isolates belonging to subtrees containing isolates from only those facilities
mixed_facility_subtree_1 <- function(subtrees, facility)
{
  
  #GET FACILITY PAIRS
  facils = sort(unique(facility));
  
  facil_pair_counts = matrix(0, ncol = length(facils), nrow = length(facils),
                             dimnames = list(facils, facils))
  
  facil_pairs = rbind(t(combn(facils,2)), t(combn(facils,2))[,2:1]);
  
  for(p in 1:nrow(facil_pairs))
  {
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVE AT LEAST ONE LABEL FOR THE PAIR OF FACILITIES OF INTEREST
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(facility[intersect(st$tip.label, names(facility))])) == 2 &&
                            facil_pairs[p,1] %in% facility[intersect(st$tip.label, names(facility))] &&
                            facil_pairs[p,2] %in% facility[intersect(st$tip.label, names(facility))])
                         {
                           
                           intersect(names(facility), st$tip.label)
                           
                         }
                       });
    
    
    facil_pair_counts[facil_pairs[p,1], facil_pairs[p,2]] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  #CALCULATE SUB-TREE SIZE FOR SINGLE FACILITY
  for(f in facils)
  {
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVING LABELS ONLY FOR SINGLE FACILITY
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(facility[intersect(st$tip.label, names(facility))])) == 1 &&
                            sum(facility[intersect(st$tip.label, names(facility))] == f) > 0 )
                         {
                           
                           intersect(names(facility), st$tip.label)
                           
                         }
                       });
    
    
    facil_pair_counts[f, f] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  return(facil_pair_counts)
  
}#end mixed_facility_subtree_1


#mixed_facility_subtree_1_rand - Takes as input: 1) subtrees          - A list produced by subtrees, that includes all isolates of interest
#                                                2) facility           - A named vactor of facility labeled
#
#For each pair of facilities returns the number of isolates belonging to subtrees containing isolates from only those facilities
#	- Randomization performed by fixing isolates from facility 1 and counting number of isolates in mixed subtrees
mixed_facility_subtree_1_rand <- function(subtrees, facility)
{
  
  #GET FACILITY PAIRS
  facils = sort(unique(facility));
  
  facil_pair_counts = matrix(0, ncol = length(facils), nrow = length(facils),
                             dimnames = list(facils, facils))
  
  facil_pairs = rbind(t(combn(facils,2)), t(combn(facils,2))[,2:1]);
  
  for(p in 1:nrow(facil_pairs))
  {
    
    #RANDOMIZE FACILITIES KEEPING ONE FACILITIES ISOLATES IN THE SAME PLACE ON THE TREE
    #THIS IS DESIGNED TO PREVENT SINGLE FACILITY CLUSTERING FROM PROVIDING DOMINANT SIGNAL
    rand_facility = facility;
    rand_facility_subset_names = names(rand_facility)[rand_facility!= facil_pairs[p,1]];
    
    names(rand_facility)[rand_facility!= facil_pairs[p,1]] = rand_facility_subset_names[sample(1:length(rand_facility_subset_names))]
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVE AT LEAST ONE LABEL FOR THE PAIR OF FACILITIES OF INTEREST
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(rand_facility[intersect(st$tip.label, names(rand_facility))])) == 2 &&
                            facil_pairs[p,1] %in% rand_facility[intersect(st$tip.label, names(rand_facility))] &&
                            facil_pairs[p,2] %in% rand_facility[intersect(st$tip.label, names(rand_facility))])
                         {
                           
                           intersect(names(rand_facility), st$tip.label)
                           
                         }
                       });
    
    
    facil_pair_counts[facil_pairs[p,1], facil_pairs[p,2]] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  #CALCULATE SUB-TREE SIZE FOR SINGLE FACILITY
  for(f in facils)
  {
    
    #RANDOMIZE FACILITIES ON TREE
    rand_facility = facility;
    names(rand_facility) = names(rand_facility)[sample(1:length(rand_facility))]
    
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVING LABELS ONLY FOR SINGLE FACILITY
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(rand_facility[intersect(st$tip.label, names(rand_facility))])) == 1 &&
                            sum(rand_facility[intersect(st$tip.label, names(rand_facility))] == f) > 0 )
                         {
                           
                           intersect(names(rand_facility), st$tip.label)
                           
                         }
                       });
    
    
    facil_pair_counts[f, f] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  return(facil_pair_counts)
  
}#end mixed_facility_subtree_1_rand


#FACILITY CLUSTERING METHOD 2
#mixed_facility_subtree_2 - Takes as input: 1) subtrees          - A list produced by subtrees, that includes all isolates of interest
#                                           2) facility          - A named vactor of facility labeled
#
#For each pair of facilities returns the number of isolates from the first first facility belonging to subtrees containing isolates from only those facilities
mixed_facility_subtree_2 <- function(subtrees, facility)
{
  
  #GET FACILITY PAIRS
  facils = sort(unique(facility));
  
  facil_pair_counts = matrix(0, ncol = length(facils), nrow = length(facils),
                             dimnames = list(facils, facils))
  
  facil_pairs = rbind(t(combn(facils,2)), t(combn(facils,2))[,2:1]);
  
  for(p in 1:nrow(facil_pairs))
  {
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVE AT LEAST ONE LABEL FOR THE PAIR OF FACILITIES OF INTEREST
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(facility[intersect(st$tip.label, names(facility))])) == 2 &&
                            facil_pairs[p,1] %in% facility[intersect(st$tip.label, names(facility))] &&
                            facil_pairs[p,2] %in% facility[intersect(st$tip.label, names(facility))])
                         {
                           
                           return(st$tip.label[facility[st$tip.label] == facil_pairs[p,1]])
                           
                         }
                       });
    
    
    facil_pair_counts[facil_pairs[p,1], facil_pairs[p,2]] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  #CALCULATE SUB-TREE SIZE FOR SINGLE FACILITY
  for(f in facils)
  {
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVING LABELS ONLY FOR SINGLE FACILITY
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(facility[intersect(st$tip.label, names(facility))])) == 1 &&
                            sum(facility[intersect(st$tip.label, names(facility))] == f) > 0 )
                         {
                           
                           intersect(names(facility), st$tip.label)
                           
                         }
                       });
    
    
    facil_pair_counts[f, f] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  return(facil_pair_counts)
  
}#end mixed_facility_subtree_2


#mixed_facility_subtree_2_rand - Takes as input: 1) subtrees          - A list produced by subtrees, that includes all isolates of interest
#                                                2) facility           - A named vactor of facility labeled
#
#For each pair of facilities returns the number of isolates from the first first facility belonging to subtrees containing isolates from only those facilities
#	- Randomization performed by fixing isolates from facility 2 and counting number of isolates from facility 1 in mixed subtrees
mixed_facility_subtree_2_rand <- function(subtrees, facility)
{
  
  #GET FACILITY PAIRS
  facils = sort(unique(facility));
  
  facil_pair_counts = matrix(0, ncol = length(facils), nrow = length(facils),
                             dimnames = list(facils, facils))
  
  facil_pairs = rbind(t(combn(facils,2)), t(combn(facils,2))[,2:1]);
  
  for(p in 1:nrow(facil_pairs))
  {
    
    #RANDOMIZE FACILITIES KEEPING ONE FACILITIES ISOLATES IN THE SAME PLACE ON THE TREE
    #THIS IS DESIGNED TO PREVENT SINGLE FACILITY CLUSTERING FROM PROVIDING DOMINANT SIGNAL
    rand_facility = facility;
    rand_facility_subset_names = names(rand_facility)[rand_facility!= facil_pairs[p,2]];
    
    names(rand_facility)[rand_facility!= facil_pairs[p,2]] = rand_facility_subset_names[sample(1:length(rand_facility_subset_names))]
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVE AT LEAST ONE LABEL FOR THE PAIR OF FACILITIES OF INTEREST
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(rand_facility[intersect(st$tip.label, names(rand_facility))])) == 2 &&
                            facil_pairs[p,1] %in% rand_facility[intersect(st$tip.label, names(rand_facility))] &&
                            facil_pairs[p,2] %in% rand_facility[intersect(st$tip.label, names(rand_facility))])
                         {
                           
                           st$tip.label[facility[st$tip.label] == facil_pairs[p,1]]
                           
                         }
                       });
    
    
    facil_pair_counts[facil_pairs[p,1], facil_pairs[p,2]] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  #CALCULATE SUB-TREE SIZE FOR SINGLE FACILITY
  for(f in facils)
  {
    
    #RANDOMIZE FACILITIES ON TREE
    rand_facility = facility;
    names(rand_facility) = names(rand_facility)[sample(1:length(rand_facility))]
    
    
    #DETERMINE THE LARGEST CLUSTER WHICH EACH ISOLATE BELONGS TO. CLUSTERS ARE DEFINED AS:
    # 1) HAVING LABELS ONLY FOR SINGLE FACILITY
    pure_tips = sapply(subtrees,
                       FUN = function(st){
                         if(length(unique(rand_facility[intersect(st$tip.label, names(rand_facility))])) == 1 &&
                            sum(rand_facility[intersect(st$tip.label, names(rand_facility))] == f) > 0 )
                         {
                           
                           intersect(names(rand_facility), st$tip.label)
                           
                         }
                       });
    
    
    facil_pair_counts[f, f] = length(unique(unlist(pure_tips)))
    
  }#end for
  
  
  return(facil_pair_counts)
  
}#end mixed_facility_subtree_2_rand


#FACILITY CLUSTERING METHOD 3
#mixed_facility_subtree_3 - Takes as input: 1) tree		 - A tree containing all isolates of interest
#					    2) subtrees          - A list produced by subtrees, that includes all isolates of interest
#                                           3) facilty          - A named vactor of facility labeled
#
#For each pair of facilities returns the number of subtrees containing only isolates from that pair of facilities
mixed_facility_subtree_3 <- function(tree, subtrees, facility)
{
  
  #GET FACILITY PAIRS
  facils = sort(unique(facility));
  
  facil_pair_counts = matrix(0, ncol = length(facils), nrow = length(facils),
                             dimnames = list(facils, facils))
  
  facil_pairs = rbind(t(combn(facils,2)), t(combn(facils,2))[,2:1]);
  
  for(p in 1:nrow(facil_pairs))
  {
    
    #DETERMINE THE NUMBER OF SUBTREES CONTAINING ISOLATES FROM BOTH FACILITY PAIRS
    pure_sts = sapply(subtrees,
                      FUN = function(st){
                        if(length(unique(facility[intersect(st$tip.label, names(facility))])) == 2 &&
                           facil_pairs[p,1] %in% facility[intersect(st$tip.label, names(facility))] &&
                           facil_pairs[p,2] %in% facility[intersect(st$tip.label, names(facility))])
                        {
                          
                          intersect(names(facility), st$tip.label)
                          
                        }
                      });
    
    #GET THE MAXIMAL SUBTREE THAT EACH ISOLATE BELONGS TO
    max_st = sapply(names(facility)[facility %in% facil_pairs[p,]],  
                    FUN = function(x)
                    {
                      #GET SIZE OF SUBTREES THAT ISOLATE BELONGS TO
                      st_length = unlist(lapply(pure_sts, 
                                                FUN = function(st_tips)
                                                {
                                                  
                                                  if(x %in% st_tips)
                                                  {
                                                    length(st_tips)
                                                  }else{
                                                    
                                                    return(0);
                                                    
                                                  }
                                                  
                                                }))
                      
                      #RETURN INDEX OF LARGEST SUBTREE
                      if(max(st_length) != 0)
                      {
                        
                        which.max(st_length);
                        
                      }else{
                        
                        return(0);
                        
                      }
                      
                    })
    
    facil_pair_counts[facil_pairs[p,1], facil_pairs[p,2]] = length(unique(max_st[max_st > 0]))
    
  }#end for
  
  return(facil_pair_counts)
  
}#end mixed_facility_subtree_3


#FACILITY CLUSTERING METHOD 4
#mixed_facility_subtree_4 - Takes as input: 1) tree              - A tree containing all isolates of interest
#					    2) subtrees          - A list produced by subtrees, that includes all isolates of interest
#                                           3) facility          - A named vactor of facility labeled
#
#For each pair of facilities returns the maximal number of subtrees (considers embedded subtrees) containing only isolates from that pair of facilities
mixed_facility_subtree_4 <- function(tree, subtrees, facility)
{
  
  #GET FACILITY PAIRS
  facils = sort(unique(facility));
  
  facil_pair_counts = matrix(0, ncol = length(facils), nrow = length(facils),
                             dimnames = list(facils, facils))
  
  facil_pairs = rbind(t(combn(facils,2)), t(combn(facils,2))[,2:1]);
  
  #COLLAPSE SUBTREES WITH ONLY ONE FACILITY PRESENT
  #tree = collapse_pure_subtrees(tree,facility)
  if(!is.binary(tree)){
    print(paste(length(table(tree$edge[,1])[table(tree$edge[,1]) != 2]), 'out of', 
                length(table(tree$edge[,1])), 'internal nodes are unresolved. Unresolved nodes will be randomly resolved.'))
    # randomly resolve polytomies
    tree = multi2di(tree,tol=0)
  }
  subtrees = subtrees(tree)
  
  for(p in 1:nrow(facil_pairs))
  {
    
    #DETERMINE THE NUMBER OF SUBTREES CONTAINING ISOLATES FROM BOTH FACILITY PAIRS
    pure_sts = sapply(subtrees,
                      FUN = function(st){
                        if(length(unique(facility[intersect(st$tip.label, names(facility))])) == 2 &&
                           facil_pairs[p,1] %in% facility[intersect(st$tip.label, names(facility))] &&
                           facil_pairs[p,2] %in% facility[intersect(st$tip.label, names(facility))])
                        {
                          
                          intersect(names(facility), st$tip.label)
                          
                        #} else {
                        #  return(0)
                        }
                      });
    
    #### WORKING HERE ####
    #GET THE MAXIMAL SUBTREE THAT EACH ISOLATE BELONGS TO
    min_st = sapply(names(facility)[facility %in% facil_pairs[p,]],  
                    FUN = function(x)
                    {
                      
                      #GET SIZE OF SUBTREES THAT ISOLATE BELONGS TO
                      st_length = unlist(lapply((pure_sts[lapply(pure_sts,length) == 2]), 
                                                FUN = function(st_tips)
                                                {
                                                  
                                                  if(x %in% st_tips)
                                                  {
                                                    length(st_tips)
                                                  }else{
                                                    
                                                    return(Inf);
                                                    
                                                  }
                                                  
                                                }))
                      

                      #RETURN INDEX OF SMALLEST SUBTREE
                      if(suppressWarnings(min(st_length)) != Inf)
                      {
                        
                        which.min(st_length);
                        
                      }else{
                        
                        return(0);
                        
                      }
                      
                    })
    
    facil_pair_counts[facil_pairs[p,1], facil_pairs[p,2]] = length(unique(min_st[min_st > 0]))
    
  }#end for
  
  return(facil_pair_counts)
  
}#end mixed_facility_subtree_4


#pure_subtree_rand - Takes as input: 1) tree       - A tree containing all isolates of interest
#                                    2) facilty    - A named vactor of facility labeled

#Creates randomized tree where pure subtrees are randomly swapped for one another
pure_subtree_rand <- function(tree, facility)
{
  
  #GET SUBTREES THAT ONLY CONTAIN ISOLATES FROM A SINGLE FACILITY
  sts = subtrees(tree);
  st_length = sapply(sts, FUN = function(st){length(st$tip.label)})
  
  pure_st_i = sapply(sts, FUN = function(st)
  {
    
    return(length(unique(facility[st$tip.label])) == 1);
    
  })
  if (sum(pure_st_i) > 0){
    
    #GET THE LARGEST PURE SUBTREE THAT EACH ISOLATE BELONGS TO
    max_st = sapply(names(facility), 
                    FUN = function(x)
                    {
                      
                      #DETERMINE WHICH SUBTREES AN ISOLATE BELONGS TO
                      st_i = sapply(sts[pure_st_i], FUN = function(st){x %in% st$tip.label})
                      
                      
                      #GET THE LARGEST SUBTREE THAT AN ISOLATE BELONGS
                      if(sum(st_i) > 0)
                      {
                        
                        st_max_i = which.max(st_length[which(pure_st_i)[st_i]]);#not getting the subtree?
                        which(pure_st_i)[which(st_i)[st_max_i]]
                        
                        
                      }else{
                        
                        return(0);
                        
                      }
                      
                    })
    
    
    #GET THE ROOT NODE FOR EACH PURE SUBTREE
    root_st = sapply(unique(max_st[max_st != 0]),
                     FUN = function(st)
                     {
                       
                       return(getMRCA(tree, sts[[st]]$tip.label))
                       
                     })
    
    #FROM EACH PURE SUBTREE, KEEP THE FIRST TIP AND REMOVE THE REST
    subtree_tips = sapply(root_st, FUN = function(x){tips(tree, node = x)})
    if (is.list(subtree_tips) == FALSE){subtree_tips = list(subtree_tips)}
    collapsed_tree = drop.tip(tree, unlist(lapply(subtree_tips, FUN = function(x){x[2:length(x)]})))
    
    #RANDOMLY SWAP ALL TIPS WITHOUT FIXING ANY
    rand_tree = collapsed_tree;
    root_st = rand_tree$tip.label
    
    root_st_perm = sample(root_st, length(root_st));
    rand_tree$tip.label = root_st_perm} else{
      
      rand_tree = tree
      root_st = rand_tree$tip.label
      
      root_st_perm = sample(root_st, length(root_st));
      rand_tree$tip.label = root_st_perm
      
    }
  
  return(rand_tree)
  
}#end pure_subtree_rand

#collapse_pure_subtrees - Takes as input: 1) tree       - A tree containing all isolates of interest
#                               2) facilty    - A named vactor of facility labeled

#Creates tree where pure subtrees are collapsed
collapse_pure_subtrees <- function(tree, facility)
{
  
  #GET SUBTREES THAT ONLY CONTAIN ISOLATES FROM A SINGLE FACILITY
  sts = subtrees(tree);
  st_length = sapply(sts, FUN = function(st){length(st$tip.label)})
  
  pure_st_i = sapply(sts, FUN = function(st)
  {
    
    return(length(unique(facility[st$tip.label])) == 1);
    
  })
  if (sum(pure_st_i) > 0){
    
    #GET THE LARGEST PURE SUBTREE THAT EACH ISOLATE BELONGS TO
    max_st = sapply(names(facility), 
                    FUN = function(x)
                    {
                      
                      #DETERMINE WHICH SUBTREES AN ISOLATE BELONGS TO
                      st_i = sapply(sts[pure_st_i], FUN = function(st){x %in% st$tip.label})
                      
                      
                      #GET THE LARGEST SUBTREE THAT AN ISOLATE BELONGS
                      if(sum(st_i) > 0)
                      {
                        
                        st_max_i = which.max(st_length[which(pure_st_i)[st_i]]);#not getting the subtree?
                        which(pure_st_i)[which(st_i)[st_max_i]]
                        
                        
                      }else{
                        
                        return(0);
                        
                      }
                      
                    })
    
    
    #GET THE ROOT NODE FOR EACH PURE SUBTREE
    root_st = sapply(unique(max_st[max_st != 0]),
                     FUN = function(st)
                     {
                       
                       return(getMRCA(tree, sts[[st]]$tip.label))
                       
                     })
    #FROM EACH PURE SUBTREE, KEEP THE FIRST TIP AND REMOVE THE REST
    subtree_tips = sapply(root_st, FUN = function(x){list(tips(tree, node = x))})
    #if (is.list(subtree_tips) == FALSE){subtree_tips = list(subtree_tips)}
    collapsed_tree = drop.tip(tree, unlist(lapply(subtree_tips, FUN = function(x){x[2:length(x)]})))
    
  }
  
   if(exists('collapsed_tree')){
     
     return(collapsed_tree)
     
   }else{
    
    return(tree)
  }
  
}#end collapse_pure_subtrees

#get_pure_subtrees_v2 - Takes as input: 1) tree       - A tree containing all isolates of interest
#                               2) facilty    - A named vactor of facility labeled

#Returns all pure subtrees
get_pure_subtrees_v2 <- function(tree, facility)
{
  
  #GET SUBTREES THAT ONLY CONTAIN ISOLATES FROM A SINGLE FACILITY
  sts = subtrees(tree);
  st_length = sapply(sts, FUN = function(st){length(st$tip.label)})
  
  pure_st_i = sapply(sts, FUN = function(st)
  {
    
    return(length(unique(facility[st$tip.label])) == 1);
    
  })
  
  if (sum(pure_st_i) > 0){
    
    #GET THE LARGEST PURE SUBTREE THAT EACH ISOLATE BELONGS TO
    max_st = sapply(names(facility), 
                    FUN = function(x)
                    {
                      
                      #DETERMINE WHICH SUBTREES AN ISOLATE BELONGS TO
                      st_i = sapply(sts[pure_st_i], FUN = function(st){x %in% st$tip.label})
                      
                      
                      #GET THE LARGEST SUBTREE THAT AN ISOLATE BELONGS
                      if(sum(st_i) > 0)
                      {
                        
                        st_max_i = which.max(st_length[which(pure_st_i)[st_i]]);#not getting the subtree?
                        which(pure_st_i)[which(st_i)[st_max_i]]
                        
                        
                      }else{
                        
                        return(0);
                        
                      }
                      
                    })
    
    #GET THE ROOT NODE FOR EACH PURE SUBTREE
    root_st = sapply(unique(max_st[max_st != 0]),
                     FUN = function(st)
                     {
                       
                       return(getMRCA(tree, sts[[st]]$tip.label))
                       
                     })
    
    pure_subtrees = lapply(root_st, function(x) extract.clade(tree, x))
    
    return(pure_subtrees)
    
    
  }
  
  
}#end get_pure_subtree_v2


