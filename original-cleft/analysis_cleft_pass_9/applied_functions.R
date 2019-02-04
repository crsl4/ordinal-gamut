library(data.table)
library(SKAT)
library(VGAM)
library(CompQuadForm)

##------------------------------------------------------------------------------
## function to generate variant set id file text
##------------------------------------------------------------------------------
write_gene_variants <- function(g, min_rv_per_gene=4) {
  tmp <- gdat[gdat$gene%in%g,]
  vids <- bim$vid[which(bim$marker%in%tmp$marker)]
  if(length(vids)>=min_rv_per_gene) {
    txt <- paste0(g,"\n",paste0(vids,"\n",collapse = ""),"END\n\n")
  } else {
    txt <- ""
  }
  return(txt)
}


##------------------------------------------------------------------------------
## function to write by_gene terminal call text
##------------------------------------------------------------------------------
by_gene_plink <- function(gene, dir) {
  txt <- paste0("cd ",dir,"by_gene\n\n",
                "mkdir ",gene,"\n",
                "cd ",gene,"\n\n",
                "../../plink --bfile ../../first_call/filtered --set ../../good_genes.set --gene ",gene," --recode A --out setdat\n\n")
  return(txt)
}

##------------------------------------------------------------------------------
## function to generate ordinal residuals (PO)
##------------------------------------------------------------------------------
getPOresiduals <- function(ordinal_vector, covariates_mat = NULL) {
  ordinal_vector <- as.ordered(ordinal_vector)
  if (is.null(covariates_mat)) {
    null_model <- vglm(ordinal_vector ~ 1,
                       family = cumulative(parallel = TRUE))    
  } else {
    null_model <- vglm(ordinal_vector ~ covariates_mat,
                       family = cumulative(parallel = TRUE))
  }
  return(residuals(null_model, "response"))
}


##------------------------------------------------------------------------------
## linear kernel gamut function
##------------------------------------------------------------------------------
linear_kernel_gamut <- function(X){
  out = vector("list", 2)
  names(out) = c("Kc", "ev_Kc")
  
  Kc = t(X) %*% X 
  ev_Kc = eigen(Kc, symmetric = T, only.values = T)$values  
  
  out$Kc = X %*% t(X)
  out$ev_Kc = ev_Kc[ev_Kc > 1e-08]
  return(out)
}

##------------------------------------------------------------------------------
## similarity projection kernel function
##------------------------------------------------------------------------------
proj_GAMuT_pheno <- function(X){
  out = vector("list", 2)
  names(out) = c("Kc", "ev_Kc")
  out$Kc = X%*%solve(t(X)%*%X,tol = 1e-17)%*%t(X) # Projection matrix
  out$ev_Kc = rep(1, ncol(X)) # Find eigenvalues
  return(out)	
}


##------------------------------------------------------------------------------
## gamut test function
##------------------------------------------------------------------------------
TestGAMuT <- function(Yc, lambda_Y, Xc, lambda_X) {
  
  ## test statistic:
  m = nrow(Yc) # number of subjects in study
  GAMuT = (1/m) * sum(sum(t(Yc) * Xc))  
  
  ## populate vector of all pairwise combination of eigenvalues
  ## from the phenotype and genotype similarity matrices:
  Z <- (as.matrix(lambda_Y)) %*% t(as.matrix(lambda_X))
  Zsort <- sort(Z, decreasing = T)
  
  ## derive p-value of GAMuT statistic:
  scoredavies = GAMuT*m^2
  results_score <- davies(scoredavies, Zsort)
  davies_pvalue <- (results_score$Qq)
  
  return(davies_pvalue)
}

##------------------------------------------------------------------------------
## perform burden, skat, and gamut gene-level tests
##------------------------------------------------------------------------------
do_gene_test <- function(gene, cov, dir, pheno_mat="linear kernel") {
  
  file <- paste0(dir,"by_gene/",gene,"/setdat.raw")
  raw <- mean_impute(file)
  
  ## check order correct
  if(!sum(raw$FID==cov$FID)==nrow(raw)) {
    stop("FIDs of cov file and raw file not in same order")
  }
  if(!sum(raw$IID==cov$IID)==nrow(raw)) {
    stop("IIDs of cov file and raw file not in same order")
  }
  
  FH <- cov$FH ## assumes FH is 0/1 variable and first covariates listed in covar file
  COVAR <- as.matrix(cov[,-1:-3])
  PHENO <- (raw$PHENOTYPE)-1
  GENO <- as.matrix(raw[,-1:-6])
  
  ## 0 = control, 1 = case w/o family history, 2 = case w/ history
  ordinal_score <- PHENO + PHENO*(FH) # family history
  
  # remove monomorphic SNPs
  GENO <- t(GENO) # rows are SNPs and columns are subjects
  GENO <- GENO[rowSums(GENO) != 0, ]
  
  # weighting (beta weights via original SKAT)
  MAF_estimate = rowMeans(GENO)/2
  beta_weight = dbeta(MAF_estimate, 1, 25) / dbeta(0, 1, 25) 
  
  # burden
  geno_burden <- colSums(beta_weight %*% GENO) # weighted sum
  model <- glm(PHENO ~ geno_burden + COVAR, family = binomial) 
  burden <- summary(model)$coefficients[2, 4]
  
  # skat
  skat_obj <- SKAT_Null_Model(PHENO ~ COVAR, out_type = "D") 
  skat <- SKAT(Z = t(GENO), obj = skat_obj, is_check_genotype = F, weights = beta_weight)$p.value
  
  # gamutlike ordinal outcome
  if (length(table(ordinal_score)) == 2) {
    resids <- getPOresiduals(PHENO)
  } else {
    resids <- getPOresiduals(ordinal_score, COVAR)
  }
  
  if(pheno_mat=="projection matrix"){
    Y <- proj_GAMuT_pheno(resids)
  } else if(pheno_mat=="linear kernel"){
    Y <- linear_kernel_gamut(resids)
  } else{
    stop("pheno_mat must be one of {'projection matrix', 'linear kernel'}")
  }
  
  G0 = t(GENO) %*% diag(beta_weight)               # weighted rare variants
  G = as.matrix(scale(G0, center = T, scale = F))  # centered genotype matrix
  linear_weighted_geno <- linear_kernel_gamut(G) 
  
  gamut <- TestGAMuT(Yc = Y$Kc, lambda_Y = Y$ev_Kc,
                     Xc = linear_weighted_geno$Kc, 
                     lambda_X = linear_weighted_geno$ev_Kc)
  
  p <- c(burden, skat, gamut)
  p <- paste(gene, paste0(p,collapse=" "))
  write(p, paste0(dir,"by_gene/gene_test_results.txt"), append = T)
  print(p)
  
  return(p)
}



##------------------------------------------------------------------------------
## function to impute missing genotypes by probability
##------------------------------------------------------------------------------
prob_impute <- function(geno_file) {

  #read genotype file and remove the first 6 columns
  data_full <- read.table(geno_file, header = T, stringsAsFactors = F)
  data <- data_full[,-c(1:6) ]
  
  if(sum(is.na(data))==0) {
    data <- cbind(data_full[, 1:6], data)
  } else {
    #calculate probabilities of the three genotypes for each column (SNP)
    n <- nrow(data)
    p_0 <- apply(data, 2, function(x){sum(x == 0, na.rm = T)/(n - sum(is.na(x)))})
    p_1 <- apply(data, 2, function(x){sum(x == 1, na.rm = T)/(n - sum(is.na(x)))})
    p_2 <- apply(data, 2, function(x){sum(x == 2, na.rm = T)/(n - sum(is.na(x)))})
    p <- data.frame(p_0, p_1, p_2)
    
    #make sure the probalilities add up to one so that it won't be a problem for the sampling function
    p <- t(apply(p, 1, function(x){x/sum(x)}))
    
    #make a table for indices of missing genotypes
    NA_indices <- which(is.na(data), arr.ind = T)
    
    #replace missing genotypes by sampling from (0,1,2) based on probabilities given in table p
    for (i in 1:nrow(NA_indices)) {
      x <- NA_indices[i, ]
      data[x[1], x[2]] <- sample(c(0:2), 1, replace = T, prob = p[x[2], ])
    }
    
    #bind the first 6 columns to the imputed data
    data <- cbind(data_full[, 1:6], data)
  }
  return(data)
}


##------------------------------------------------------------------------------
## function to impute missing genotypes by mean
##------------------------------------------------------------------------------
mean_impute <- function(geno_file) {
  
  #read genotype file and remove the first 6 columns
  data_full <- read.table(geno_file, header = T, stringsAsFactors = F)
  data <- data_full[,-c(1:6) ]
  
  if(sum(is.na(data))==0) {
    data <- cbind(data_full[, 1:6], data)
  } else {
    
    #calculate average genotype for each column (SNP)
    avg <- apply(data, 2, function(x) mean(x, na.rm = T))
    
    #make a table for indices of missing genotypes
    NA_indices <- which(is.na(data), arr.ind = T)
    
    #replace missing genotypes with SNP average genotype
    for (i in 1:nrow(NA_indices)) {
      x <- NA_indices[i, ]
      data[x[1], x[2]] <- avg[x[2]]
    }
    
    #bind the first 6 columns to the imputed data
    data <- cbind(data_full[, 1:6], data)
  }
  return(data)
}



##------------------------------------------------------------------------------
## function to grab gene data
##------------------------------------------------------------------------------
get_gene_dat <- function(gene){
  rv_count <- length(which(gdat$gene%in%gene))
  snps <- paste(gdat$marker[which(gdat$gene%in%gene)],collapse="$")
  chroms <- paste(gdat$chrom[which(gdat$gene%in%gene)],collapse="$")
  unique_chroms <- unique(gdat$chrom[which(gdat$gene%in%gene)])
  if(length(unique(chroms))==1){
    loc <- substr(unique_chroms,4,nchar(unique_chroms))
  } else{
    loc <- "not unique"
  }
  return(c(rv_count, snps, chroms, loc))
}


##------------------------------------------------------------------------------
## function to write cov data
##------------------------------------------------------------------------------
write_covs <- function(i) {
  ind <- which(dat_ids%in%i)
  out <- paste0(dat[ind,c("fid","iid","Multiplex","CAUC_EV1","CAUC_EV2","CAUC_EV3","CAUC_EV4","CAUC_EV5")],collapse=" ")
  out <- paste0(out,"\n")
  return(out)
}