##------------------------------------------------------------------------------
## function to generate params file for COSI
## note: only european and african haplotypes considered
##------------------------------------------------------------------------------
writeParams <- function(irep, length, n_hap_ceu, n_hap_yri){
  print(paste0("Writing params file with length (in bp) = ", length, ", european sample size = ", 
               n_hap_ceu, ", and african sample size = ", n_hap_yri))
  text = paste0("# sample file\n",
                "# comments have #s in front of them\n",
                "# newlines don't matter.\n",
                "\n",
                "#-- options that could be uncommented\n",
                "#infinite_sites yes\n",
                "#random_seed 12345   # Specifies a particular random number seed\n",
                "#number_mutation_sites 2000  # Fixed number of mutation sites\n",
                "\n",
                "# in bp.\n")
  text = paste0(text,"length ", length,"\n",
                "\n",
                "# per bp per generation\n",
                "mutation_rate 1.5e-8\n",
                "\n",
                paste0("recomb_file model",irep,".test\n"),
                "gene_conversion_rate 4.5e-9\n",
                "\n",
                "\n",
                "# population info\n",
                "\n",
                "pop_define 1 european\n",
                "pop_define 3 african-american\n",
                "pop_define 4 asian\n",
                "pop_define 5 african\n",
                "\n",
                "#european\n",
                "pop_size 1 100000\n")
  text = paste0(text,"sample_size 1 ",n_hap_ceu,"\n",
                "\n",
                "#african american\n",
                "pop_size 3 100000\n",
                "sample_size 3 0\n",
                "\n",
                "#asian\n",
                "pop_size 4 100000\n",
                "sample_size 4 0\n",
                "\n",
                "#african\n",
                "pop_size 5 100000\n")
  text = paste0(text, "sample_size 5 ",n_hap_yri,"\n",
                "\n",
                "pop_event migration_rate \"afr->eur migration\" 5 1 0. .000032\n",
                "pop_event migration_rate \"eur->afr migration\" 1 5 0 .000032\n",
                "pop_event migration_rate \"afr->as migration\" 5 4 0. .000008\n",
                "pop_event migration_rate \"as->afr migration\" 4 5 0 .000008\n",
                "pop_event admix \"african american admix\" 3 1 5. .2\n",
                "pop_event split \"african to aa\" 5 3 7.0\n",
                "\n",
                "pop_event change_size \"agriculture - african\" 5 200 24000\n",
                "pop_event change_size \"agriculture - european\" 1 350 7700\n",
                "pop_event change_size \"agriculture - asian\" 4 400 7700\n",
                "pop_event bottleneck \"african bottleneck\" 5 1997 .008\n",
                "pop_event bottleneck \"asian bottleneck\" 4 1998 .067\n",
                "pop_event bottleneck \"european bottleneck\" 1 1999 .02\n",
                "\n",
                "pop_event split \"asian and european split\" 1 4 2000\n",
                "pop_event migration_rate \"afr->eur migration\" 5 1 1996 0\n",
                "pop_event migration_rate \"eur->afr migration\" 1 5 1995 0\n",
                "pop_event migration_rate \"afr->as migration\" 5 4 1994 0\n",
                "pop_event migration_rate \"as->afr migration\" 4 5 1993 0\n",
                "\n",
                "pop_event bottleneck \"OoA bottleneck\" 1 3499 .085\n",
                "pop_event split \"out of Africa\" 5 1 3500\n",
                "\n",
                "pop_event change_size \"african pop size\" 5 17000 12500")
  write(text,file=paste0("params",irep))
}


##------------------------------------------------------------------------------
## function to generate recParams file for COSI
## note: based on COSI best fit model
##------------------------------------------------------------------------------
writeRecParams <- function(irep){
  text = paste0("model 1\n",
                "#random_seed 56789\n",
                paste0("outfile model",irep,".test\n"),
                "distribution autosomes_deCODE.distr\n",
                "distance_shape 0.35\n",
                "local_shape 0.3\n",
                "intensity_shape 0.3\n",
                "space 8500\n",
                "bkgd 0.12\n",
                "#  The following option introduces variation into the base\n",
                "#   recombination rate within the region; specifically, the\n", 
                "#   rate changes every 100 kb. (The regional rate, picked from\n",
                "#   the deCODE map, remains the same.)\n",
                "#local_size 100000\n"
                )
  write(text,file=paste0("recParams",irep))
}



##------------------------------------------------------------------------------
## function to simulate one proband and their 1st degree family
##------------------------------------------------------------------------------
simulate_one_individual <- function(seed,p_yri,Nsibsize,include_parents,randomize_nsibs,
                                    cconfounder_family,hap1,hap5=NA,intercept,beta_cconfounder,
                                    beta_yri,which_rare_causal,rarebeta,lambda,minaff) {
    set.seed(seed)

    ## covariates
    yri <- rbinom(n = 1, size = 1, prob = p_yri)
    cconfounder <- rnorm(n = 1, mean = 0, sd = 1)

    ## random sibling generation
    Nsibsize_eff <- Nsibsize + include_parents*2
    if(randomize_nsibs == 1) {
        Nsibsize <- rpois(n = 1, lambda = 2)
        Nsibsize_eff <- Nsibsize + 2*include_parents
    }

    if(cconfounder_family == T) {
        cconfounder <- c(cconfounder, rnorm(n = Nsibsize_eff, mean = cconfounder, sd = 0.1))
    }

    ## sample from parental haplotypes
    if (yri == F) {
        hap <- hap1
    } else {
        hap <- hap5
    }
    
    hap_mom_idx <- sample(seq(1, nrow(hap)), size = 2, replace = T)
    hap_dad_idx <- sample(seq(1, nrow(hap)), size = 2, replace = T)

    ## bernoulli re: which allele (#1 or #2) inherited from mom or dad
    ## perform for each of the siblings
    trans_mom <- rbinom(Nsibsize + 1, 1, 0.5) + 1
    trans_dad <- rbinom(Nsibsize + 1, 1, 0.5) + 1
    hap_sibs <- hap[ hap_mom_idx[trans_mom], ] + hap[ hap_dad_idx[trans_dad], ]

    if (include_parents == 1) {
        hap_sibs <- rbind(hap_sibs, colSums(hap[hap_mom_idx, ]), colSums(hap[hap_dad_idx, ]))
    }
    
    ## probability of disease
    if (is.null(dim(hap_sibs))) { 
        pdis <- intercept + cconfounder*beta_cconfounder + yri*beta_yri + t(hap_sibs[which_rare_causal] %*% rarebeta)
    } else { 
        pdis <- intercept + cconfounder*beta_cconfounder + yri*beta_yri + t(hap_sibs[, which_rare_causal] %*% rarebeta)
    }
    
    pdis <- exp(pdis) / ( 1 + exp(pdis) )
    
    ## whether affected by disease
    affsib <- sapply(pdis, function(x) rbinom(n = 1, size = 1, prob = x))
    
    if (affsib[1] == 1) {
        pdissib <- pdis*lambda # NOTE: assuming same lambda btwn parents + siblings
        pdissib <- 1*(pdissib > 1) + pdissib*(pdissib < 1)
        affsib[-1] <- sapply(pdissib, function(x) rbinom(n = 1, size = 1, prob = x))[-1]
    } else {
      pdissib <- rep(NA,Nsibsize+1) # for flagged info #############
    }

    ## output
    disease <- affsib[1]
    num_family_affected <- sum(affsib[-1])

    ## 0 = control, 1 = case w/o family history, 2 = case w/ history
    ordinal_score <- disease + disease*(num_family_affected >= minaff) # family history
    
    if (is.null(dim(hap_sibs)) == F) {
        hap_sibs <- hap_sibs[1, ]
    }
    
    # return (6 + number of variants) variables
    return(c(ordinal_score, disease, num_family_affected, yri, cconfounder[1], seed, hap_sibs))
}


##------------------------------------------------------------------------------
## function to generate ordinal residuals (proportional odds)
##------------------------------------------------------------------------------
getPOresiduals <- function(ordinal_vector, covariates_vector = NULL) {
  ordinal_vector <- as.ordered(ordinal_vector)
  if (is.null(covariates_vector)) {
    null_model <- vglm(ordinal_vector ~ 1,
                       family = cumulative(parallel = TRUE))    
  } else {
    null_model <- vglm(ordinal_vector ~ covariates_vector,
                       family = cumulative(parallel = TRUE))
  }
  return(residuals(null_model, "response"))
}


##------------------------------------------------------------------------------
## function to generate ordinal residuals (continuation ratio)
##------------------------------------------------------------------------------
getCRresiduals <- function(ordinal_vector, covariates_vector) {
  ordinal_vector <- as.ordered(ordinal_vector)
  null_model <- vglm(ordinal_vector ~ covariates_vector,
                     family = cratio(parallel = TRUE))
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
## function to perform burden, skat, and gamut tests
##------------------------------------------------------------------------------
method_comparisons <- function(sim_results, residual_format) {
  sim_results <- as.matrix(sim_results)
  ordinal_score <- sim_results[1,]
  disease <- sim_results[2,]
  num_family_affected <- sim_results[3,]
  yri <- sim_results[4,]
  cconfounder <- sim_results[5,]
  geno <- sim_results[7:nrow(sim_results),]
  
  # remove monomorphic SNPs
  geno <- geno[rowSums(geno) != 0, ]
  
  # weighting (beta weights via original SKAT)
  MAF_estimate = rowMeans(geno)/2
  beta_weight = dbeta(MAF_estimate, 1, 25) / dbeta(0, 1, 25) 
  
  # burden
  geno_burden <- colSums(beta_weight %*% geno) # weighted sum
  model <- glm(disease ~ geno_burden, family = binomial)
  burden <- summary(model)$coefficients[2, 4]
  
  # skat
  skat_obj <- SKAT_Null_Model(disease ~ 1, out_type = "D")
  skat <- SKAT(Z = t(geno), obj = skat_obj, is_check_genotype = F,
               weights = beta_weight)$p.value
  
  # gamutlike, ordinal outcome
  if (length(table(ordinal_score)) == 2) {
    if (residual_format == "PO") { resids <- getPOresiduals(disease) }
    if (residual_format == "CR") { resids <- getCRresiduals(disease) }
  } else {
    if (residual_format == "PO") { resids <- getPOresiduals(ordinal_score, cbind(cconfounder)) }
    if (residual_format == "CR") { resids <- getCRresiduals(ordinal_score, cbind(cconfounder)) }
  }
  
  Y <- linear_kernel_gamut(resids)
  
  G0 = t(geno) %*% diag(beta_weight)               # weighted rare variants
  G = as.matrix(scale(G0, center = T, scale = F)) # centered genotype matrix
  linear_weighted_geno <- linear_kernel_gamut(G) 
  
  gamut <- TestGAMuT(Yc = Y$Kc, lambda_Y = Y$ev_Kc,
                     Xc = linear_weighted_geno$Kc, lambda_X = linear_weighted_geno$ev_Kc)
  
  return( c(burden, skat, gamut) )
}


##------------------------------------------------------------------------------
## function to calculate variant frequencies per group
##------------------------------------------------------------------------------
variant_statistics_freq_DT <- function(geno_table, causal_vector, genosAsRows=TRUE) {
  
  if ( length(geno_table) == 0) {
    return("-999,-999,-999") #return("none in simulated pop,,")
  }
  
  
  if ( is.null(dim(geno_table)) ) {
    return("-999,-999,-999") }
  
  # genos as rows
  if (genosAsRows == T) {
    
    # genosAsRows = T
    tmp <- c(
      # average MAF
      mean( rowMeans(geno_table) / 2 ),
      # average number of variants
      mean( rowSums(geno_table != 0) ),
      # average number of causal variants
      mean( rowSums(geno_table[causal_vector, ] != 0) ))
  }
  # genos as columns
  if (genosAsRows == F) {
    # genosAsRows = F
    tmp <- c(
      # average MAF
      mean( colMeans(geno_table) / 2 ),
      # average number of variants
      mean( colSums(geno_table != 0) ),
      # average number of causal variants
      mean( colSums(geno_table[, causal_vector] != 0) ))
  }
  
  return(paste(tmp, collapse = ", "))
}

