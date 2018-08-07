# Power calculator from:
# A phylogeny-based sampling strategy and power calculator informs genome-wide associations study design for microbial pathogens
#Farhat et al. Genome medicine 2014
# https://doi.org/10.1186/s13073-014-0101-7

########      Site Level Power Calculator Function
########      based on binomial distribution

########      Input Variables     ##############
## k     = Genome length in bp
## n     = vector of sample size (number of strain pairs)
## s     = vectors of pairwise distance between matched pairs
## np    = nominal false positive rate (default 0.05) this will be bonferroni corrected
## f     = vector of effect sizes to test (range 0-1)

#######       Output variables    ##############
## l     = threshold number of mutations in ph+ vs ph- strains to call a site significant
##         this is also the quantile of null distribution that protects the nominal fp rate
## power = true positive rate in percent
## fp    = exact false positive rate at the threshold l

powerCalculatorSite <- function(k,s,f,n,np)
{
  dimnames1<-list(NULL,c("n","s","l","f", "power","fp"))
  powermatrixsite<-matrix(nrow=1, ncol=6, dimnames=dimnames1, data=NA)
  if (missing(np)){
    np<-0.05
  } 
  for (e in f) {
    for (m in n) {
      for (r in s) {
        ts<-m*r #number of snp sites to be tested assumes that k>>> n*s and so overlap between snps is unlikely under neutral evolution
        fp<-np/ts ##nominal false positive rate bonferroni corrected
        p<-r/k
        q<-qbinom(fp,m,p, lower.tail=FALSE)  ##quantile of null that protects nominal fp rate
        efp<-pbinom(q,m,p, lower.tail=FALSE) ##actual fp rate
        tp<-pbinom(q,m,e, lower.tail=FALSE)  ##probability of true positive at the quantile given the alternative
        powermatrixsite<-rbind(powermatrixsite, c(m, r, q, e, tp*100, efp))
      }
    }
  }
  powermatrixsite<-powermatrixsite[-1,]
  return(as.data.frame(powermatrixsite))
}


########      Locus Level Power Calculator Function
########      based on poisson distribution

########      Input Variables     ##############
## k     = Genome length in bp
## n     = vector of sample size (number of strain pairs)
## s     = vectors of pairwise distance between matched pairs
## np    = nominal false positive rate (default 0.05) this will be bonferroni corrected
## f     = vector of effect sizes to test (range 0-1)
## gl    = vector of locus lengths for the genome of interest

#######       Output variables    ##############
## l     = threshold number of mutations in ph+ vs ph- strains to call a site significant
##         this is also the quantile of null distribution that protects the nominal fp rate
## power = true positive rate in percent
## fp    = exact false positive rate at the threshold l


powerCalculatorLocus <- function(k,s,f,n,gl,np)
{
  gl<-sort(gl[[1]]) #sort loci by length
  glq<-quantile(gl,probs=seq(0,1,0.01)) ##seperate loci into 100 quantiles
  dimnames1<-list(NULL,c("n","s","l","f","power","fp"))
  powermatrix<-matrix(nrow=1, ncol=6, dimnames=dimnames1, data=NA)
  if (missing(np)){
    np<-0.05
  } 
  for (e in f) {
    for (m in n) {
      for (x in s) {
        efp<-0
        tp<-0
        for (i in seq(2,length(glq),1)) {
          g<-(glq[[i]]+glq[[i-1]])/2
          p<-m*x*g/k
          tl<-k*(1-exp(-p))/g  #number of loci to be tested
          fp<-0.05/tl  ##bonferonni adjusted p for number of tests     
          q<-qpois(fp, p, lower.tail=FALSE)
          tefp<-ppois(q,p, lower.tail=FALSE)  ##actual fp rate
          ttp<-ppois(q,m*e, lower.tail=FALSE)  ##probability of true positive at the quantile given the alternative
          efp<-efp+tefp*0.01 #integrating the false positive rate
          tp<-tp+ttp*0.01  #integrating the true positive rate 
        }
        powermatrix<-rbind(powermatrix, c(m, x, q, e,tp*100, efp))  
      }
    }
  }
  powermatrix<-powermatrix[-1,]
  return(as.data.frame(powermatrix))
}

##Plotting function for calculator results
## b          vector to plot by options are f, s
## "nameby"   name of the variable vector to plot by options are f, s
## fix        vector of the variable to fix options are f, s
## "namefix"  name of the variable to fix options are f, s
## value      the numerical value to fix variable fix to eg s=0.20  ##plotter will return error if this value was not testing in the calculation
## level      locus or site
plotPowerResults <-function(result,b,nameb, fix, namefix, value, level)  
{
  pm<-result[which(fix==value),]
  b<-b[which(fix==value)]
  plot(pm$n, pm$power, main=paste(level," Level Power by ", nameb, "\n", namefix, "=", value,sep=""), xlab="sample size n-pairs", ylab="power (percent)", type="n")
  palette<-c("white", "black", "brown", "darkred", "red", "orange", "green", "blue")
  cl<-colors()
  cl<-cl[-c(153:252)]
  v<-unique(b)
  for (i in 1:length(v)) {
    pm1<-pm[b==v[i],]
    if (i>7) {
      ii<-i*5
    } else {
      ii<-i
    }
    palette<-c(palette, cl[ii])
    lines(pm1$n, pm1$power, col=palette[i+1])
  }
  legend(x=240, y=50, legend=c(paste(nameb, sep=""),v), col=palette, lty=1, cex=0.75, seg.len=0.6)
}