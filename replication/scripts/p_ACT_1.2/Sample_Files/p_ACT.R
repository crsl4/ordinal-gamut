# p_ACT.R 
# version 1.1  5/29/2007 (comments updated 7/23/2007)
# version 1.2  2/22/2008 Minor fixes made to p_ACT_seq.R.  No changes made to p_ACT.R.

version=1.2

library(MASS)
library(mvtnorm)
# If not already installed, the package 'mvtnorm' may be downloaded at
# http://cran.r-project.org/src/contrib/Descriptions/mvtnorm.html
# For > 100 tests, version 0.8-0 of mvtnorm (or higher) is required.

sink("p_ACT.out")
cat("p_ACT.R, version",version,"\n")
cat("************************************************\n")

err=0  #Parameter to keep track of error in input files

# These parameters affect precision of p-values and speed of estimates.  The initial settings request
# increased precision for values of p_ACT <.05, <.01, and <.00001 but may be altered as desired.
# Use 25000 for quick, reasonably precise p-values.  Use 25000000 for ultra-precise but slower p-values.
level1=25000
cutoff2=.05
level2=250000
cutoff3=.01
level3=2500000
cutoff4=.00001
level4=25000000

# The following commands read in user-created (flat AASCI space-delimited) datafiles 
# containing information on L association tests based on K traits and M markers.  
# Formatting details in README.txt.

# pvalues.txt should contain L rows and either 3 or 4 columns.
# If all tests are two-sided, only 3 columns are necessary:
#      Column 1: trait name;   Column 2: SNP name;   Column 3: p-value.
# If any one-sided tests were performed, a fourth column is needed:
#      Column 4: -1 or 1 for one-sided test, 2 for two-sided test
#        (-1 vs. 1 indicates direction of result; see README.txt for more information.)
pvals=read.table('pvalues.txt')
L=dim(pvals)[1]

# genotype.txt contains a row of M genotype scores (e.g. allele counts) for each individual, 
# plus a header row containing marker names, with missing values coded as NA.
genolist=pvals[1,2]
M=corg=1
if (file.access('genotype.txt')==0) {
   G=data.matrix(read.table('genotype.txt',header=TRUE))
   N=dim(G)[1]
   M=dim(G)[2]
   
   if (M>1) {
      #For covariance estimation, fill in missing values with mean genotype code
      for (i in 1:M) { G[is.na(G[,i]),i]=mean(G[,i],na.rm=TRUE) }

      # If there is no covar.txt file, p_act will be computed assuming there are no covariates.
      # If the model contains covariates, covar.txt should contain a row of covariates for each individual, plus
      # a header row containing covariate names, with missing values coded as NA.  See README.txt for more info.
      X=matrix(1,N)

      if (file.access('covar.txt')==0) {
          covar=data.matrix(read.table('covar.txt',header=TRUE))
          if (dim(covar)[1]!=N) {
             cat("Error: Different # of individuals in genotype.txt and covar.txt \n")
             err=1
          }
          if (dim(covar)[1]==N) {
             X=cbind(X,covar)
             #For covariance estimation, fill in missing values with variable mean
             for (i in 1:dim(X)[2]) { X[is.na(X[,i]),i]=mean(X[,i],na.rm=TRUE) }
          }
      }
      #Computation of variance matrix of G, conditional on X
      genolist=t(read.table('genotype.txt',header=FALSE)[1,])
      kg=order(genolist)
      G=G[,kg]
      corg=cov2cor(t(G)%*%G - t(G)%*%X %*% ginv(t(X)%*%X) %*% t(X)%*%G)
   }
}

# traits.txt contains a row of K trait values for each individual, plus a header row containing trait
# names, with missing values coded as NA.  If the model includes environmental covariates, trait
# values should be residualized on covariates (see README.txt for more details).
# If there is no traits.txt file, it will be assumed that only a single trait is considered.
traitlist=pvals[1,1]
K=cory=1
if (file.access('traits.txt')==0) {
   Y=data.matrix(read.table('traits.txt',header=TRUE))
   if (exists("N")) {
      if (dim(Y)[1]!=N) {
          cat("Error: Different # of individuals in genotype.txt and traits.txt \n")
          err=1
      }
   }
   if (!exists("N")) {N=dim(Y)[1]}
   K=dim(Y)[2]
   if (K>1) {
      #For covariance estimation, fill in missing values with trait mean
      for (i in 1:K) { Y[is.na(Y[,i]),i]=mean(Y[,i],na.rm=TRUE) }
      #Sort columns of Y into order by trait name, compute correlation matrix for Y
      traitlist=t(read.table('traits.txt',header=FALSE)[1,])
      kl=order(traitlist)
      Y=Y[,kl]
      cory=cor(Y)
   }
}

if (K*M<L)  {
    cat("Error: More p-values than trait*genotype combinations \n")
    err=1
}

if (err==0) {    #Condition on err==0

#Computation of correlation matrix between K x M tests
v=kronecker(cory,corg)

#Removal of any missing tests from correlation matrix
if (K*M > L) {
    glist=cbind(1,as.matrix(genolist))
    tlist=cbind(1,as.matrix(traitlist))
    tests_km=merge(tlist,glist,by.x=1,by.y=1)[,2:3]
    all=merge(tests_km,pvals,by.x=c(1,2),by.y=c(1,2),all=TRUE)
    v=v[!is.na(all[,3]),!is.na(all[,3])]
}

#Computation of p_ACT
minp=min(pvals[,3])
rank=order(pvals[,3])
ordered=pvals[rank,]
best=ordered[1,]

if (minp==0) {p_ACT=0}
if (minp>=.5) {p_ACT=1}

if (minp>0 & minp < .5) {
   lower=rep(qnorm(minp/2),L)
   upper=rep(qnorm(1-minp/2),L)
   if (dim(pvals)[2]>3) {
       p4=pvals[,4]
       lower=rep(-Inf,L)
       upper=rep(Inf,L)
       lower[p4==2]=qnorm(minp/2)
       lower[p4==-1]=qnorm(minp)
       upper[p4==2]=qnorm(1-minp/2)
       upper[p4==1]=qnorm(1-minp)
   }
# Use pmvnorm() function from mvtnorm package to compute multivariate normal probabilities
# For more information: 
# Genz A, Bretz F, Hothorn T (2006) mvtnorm: Multivariate normal and t distribution. R package version 0.8-0
# Genz A (1992) Numerical computation of multivariate normal probabilities.  J Comput Graph Stat 1:141-150
# Genz A (1993) Comparison of methods for the computation of multivariate normal probabilities.  Comput Sci Stat 25:400-405

   p_ACT=1-pmvnorm(lower=lower,upper=upper,sigma=v,maxpts=level1,abseps=.0000000000001)
   if (p_ACT<cutoff2) {
       p_ACT=1-pmvnorm(lower=lower,upper=upper,sigma=v,maxpts=level2,abseps=.0000000000001)
       if (p_ACT<cutoff3) {
           p_ACT=1-pmvnorm(lower=lower,upper=upper,sigma=v,maxpts=level3,abseps=.0000000000001)
           if (p_ACT<cutoff4) {
               p_ACT=1-pmvnorm(lower=lower,upper=upper,sigma=v,maxpts=level4,abseps=.0000000000001)
           }
       }
   } 
}
bestlist=paste(best[1,1],best[1,2])
cat(L,"tests performed. \n")
cat("Trait and SNP with smallest p-value: ",bestlist,"\n")
cat("Smallest p-value: ",best[,3],"\n")
cat(c("p_ACT:",format(p_ACT,digits=5)),"\n \n \n")

#Remove comment characters from following 2 lines to print covariance matrix
#cat("Covariance matrix from",format(L,justify="right",trim=TRUE),"tests","\n")
#write.table(format(v,digits=2),quote=FALSE,col.names=FALSE,row.names=FALSE)

} #End condition on err==0
sink()
