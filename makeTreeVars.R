# load libraries
library(phytools)
library(ape)

# makeTreeVars - plots tree and heatmap of variants side by side. The first column is mic.
# Input:
#     tree - phylogenetic tree. tip labels must be a subset of column names in varMat
#     varMat - variance matrix with isolates named as id_mic
#     genomicRegions - matrix of positions to iterate over. nx2 matrix with col1=start, col2=end. rownames = name of region
#     pos - positions of each variant in varMat (length(pos) == length(rownames(varMat)))
#     varNames - names for columns of matrix. same length as number of rows in varMat (ex. paste(parsed$pos, parsed$var_nucleotides))
#     matType - type of matrix (snp, kmer, etc.). default: snpmat
# returns plots for each row in genomicRegions and list of variant matrices for each row in genomicRegions (genomic locations of interest)
makeTreeVars = function(tree,varMat,genomicRegions,pos,varNames=NULL,matType="snpmat",ret=F){
  
  names(genomicRegions) = c('start','end')
  
  vars = list()
  
  # ITERATE OVER EACH GENOMIC REGION OF INTEREST
  for(i in 1:length(rownames(genomicRegions))){
    #print(i)
    print(rownames(genomicRegions)[i])
    # SUBSET MATRIX TO INCLUDE ONLY VARIANTS IN GENOMIC REGION OF INTEREST
    vars[[i]] = varMat[which(pos >= genomicRegions$start[i] & pos <= genomicRegions$end[i]),]
    print(dim(vars[[i]]))
    
    if(dim(vars[[i]])[1]!=0){
    #if(!is.null(dim(vars[[i]]))){
      colnames(vars[[i]]) = colnames(varMat)
      rownames(vars[[i]]) = varNames[which(pos >= genomicRegions$start[i] & pos <= genomicRegions$end[i])]
      vars[[i]] = subset(vars[[i]],select=c(names(vars[[i]]) %in% tree$tip.label))
      vars[[i]] = unique(vars[[i]])
      mics = log(as.numeric(gsub('^.*_','',names(vars[[i]])))) #log mics so you can see differences better
      
      # PLOT JUST VARIANTS BY INCREASING MIC
      pdf(paste('figures/', rownames(genomicRegions)[i],'_',matType,'Vars.pdf',sep=''))
      varPlot(rbind(vars[[i]],mics),title=rownames(genomicRegions)[i])
      dev.off()
      
      # PLOT TREE AND VARIANTS
        pdf(paste('figures/', rownames(genomicRegions)[i],'_',matType,'TreeVars.pdf',sep=''))
        phylo.heatmap(tree,cbind(mics,t(vars[[i]])),fsize=c(0.1,0.3,1))
        title(rownames(genomicRegions)[i])
        dev.off()
      #}
    }
  } 
  if(ret) return(vars)
}

# MODIFIED FROM THE INTERNET
# Called by makeTreeVars
# ----- Define a function for plotting a matrix ----- #
varPlot <- function(x, kmers=T, thresh=NULL, ...){
  min <- min(x, na.rm=T)
  max <- max(x,na.rm=T)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  #ColorRamp <- rgb( seq(0,1,length=256),  # Red
  #                  seq(0,1,length=256),  # Green
  #                  seq(1,0,length=256))  # Blue
  ColorRamp = heat.colors(n = 2000)[2000:1]
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # if(kmers){
  #   ColorLevels <- c(0,1)
  # } else{
  #   ColorLevels <- c(-1,0,1)
  # }
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), xlab="", col=ColorRamp,
        ylab="variants", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  if(!kmers){
    axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
         cex.axis=0.7)
  }
  if(!is.null(thresh)) abline(h=thresh)
  
  
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        #xlab="",ylab="",
        xaxt="n")#, yaxt='n')
  # if(kmers){
  #   axis(LEFT<-2, at=0:1, labels=c('0','1'), las= HORIZONTAL<-0, cex.axis=0.7)
  # } else{
  #   axis(LEFT<-2, at=-1:1, labels=c('-1','0','1'), las= HORIZONTAL<-0, cex.axis=0.7)
  # }
  
  layout(1)
}
# ----- END plot function ----- #
