
## current simulation in:
args <- commandArgs(TRUE)
if(length(args) == 0){
  irep <- 1
}else{
  irep <- as.double(args[1])
}


##------------------------------------------------------------------------------
## external package dependencies:
##------------------------------------------------------------------------------
library(data.table)
library(SKAT)
library(VGAM)
library(CompQuadForm)


##------------------------------------------------------------------------------
## R code dependencies:
##------------------------------------------------------------------------------
source("00_child_functions.r")
load("arguments.Rda")


##------------------------------------------------------------------------------
## output file(s):
##------------------------------------------------------------------------------
outfile = paste0("ordinal_simulation_",irep,".txt")


##------------------------------------------------------------------------------
## create the params and recParams files 
## assumes that all individuals are european or african for now
##------------------------------------------------------------------------------
writeParams(irep, length=seqlength, n_hap_ceu=nhap_eur, n_hap_yri=nhap_afr)
writeRecParams(irep)


##------------------------------------------------------------------------------
## simulate haplotype pool:
##------------------------------------------------------------------------------
starttime <- proc.time()
cosi.log <- system(paste0("perl run.pl ",irep), intern=TRUE)

## saving and setting up seeds:
r.seed <- sample(1:555555,1,replace=FALSE)
recomb.seed <- substr(cosi.log[1],27,nchar(cosi.log[1]))
coal.seed <- substr(cosi.log[2],18,nchar(cosi.log[2]))
set.seed(r.seed)

## read in ceu snp data from cosi
pos1 <- matrix(scan(paste0("out",irep,".pos-1"), skip = 1, quiet = T), ncol = 7, byrow = T) # ceu
nvariants <- nrow(pos1)
colnames(pos1) <- c("SNP", "CHROM", "CHROM_POS", "ALLELE1", "FREQ1", "ALLELE2", "FREQ2")

hap1 <- matrix(scan(paste0("out",irep,".hap-1"), quiet = T), ncol = 2 + nvariants, byrow = T)
hap1 <- hap1[, -c(1, 2)]
hap1 <- abs(hap1 - 2) # 1/2 to 1/0 coding

## read in yri snp data from cosi for simulations with confounding (needs to be updated)
if(!p_yri==0) {
    pos5 <- matrix(scan('out.pos-5', skip = 1, quiet = T), ncol = 7, byrow = T) # yri
    colnames(pos5) <- c("SNP", "CHROM", "CHROM_POS", "ALLELE1", "FREQ1", "ALLELE2", "FREQ2")

    hap5 <- matrix(scan('out.hap-5', quiet = T), ncol = 2 + nvariants, byrow = T)
    hap5 <- hap5[, -c(1, 2)]
    hap5 <- abs(hap5 - 2) # 1/2 to 1/0 coding

    avg_maf <- rowMeans(cbind(pos1[,5],pos5[,5]))
} else {
    avg_maf <- pos1[,5]
}

## find rare variants
rarelist <- which(avg_maf <= rare_maf)
which_rare_causal <- sample(rarelist, length(rarelist)*p_causal)
maf_causal <- avg_maf[which_rare_causal]
beta_geno <- log(effectsize)/4
rarebeta <- -beta_geno*abs(log10(maf_causal)) + rnorm(n = length(maf_causal), mean = 0, sd = 0.1) # protective variants
to_reverse <- sample(x = seq(which_rare_causal),size = length(which_rare_causal)*(p_adverse)) # DOUBLE CHECK ERROR!!!!
rarebeta[to_reverse] <- -rarebeta[to_reverse] 

## calculate intercept based on target prevalence
intercept <- log(prev_target/(1 - prev_target)) - beta_yri*p_yri - sum(maf_causal*rarebeta)


##------------------------------------------------------------------------------
## simulate proband genotype and disease data:
##------------------------------------------------------------------------------
sim_results <- {}
seed <- r.seed
n_control <- 0; n_case_noFH <- 0; n_case_FH <- 0

if(p_yri==0) {
  hap5 <- NA
}

while( n_control < n_per_group | n_case_noFH < n_per_group | n_case_FH < n_per_group ) {
  sim_results_tmp <- simulate_one_individual(seed=seed, p_yri=p_yri, Nsibsize=Nsibsize, include_parents=include_parents, 
                                             randomize_nsibs=randomize_nsibs, cconfounder_family=cconfounder_family, 
                                             hap1=hap1, hap5=hap5, intercept=intercept, beta_cconfounder=beta_cconfounder,
                                             beta_yri=beta_yri, which_rare_causal=which_rare_causal, rarebeta=rarebeta, lambda=lambda, 
                                             minaff=minaff)
  
  if(sim_results_tmp[1]==0 & n_control < n_per_group) {
    sim_results <- cbind(sim_results, sim_results_tmp)
    n_control <- n_control+1
    #kept <- c(kept, 1) # for flagged info

  } else if(sim_results_tmp[1]==1 & n_case_noFH < n_per_group) {
    sim_results <- cbind(sim_results, sim_results_tmp)
    n_case_noFH <- n_case_noFH+1
    #kept <- c(kept, 1) # for flagged info

  } else if(sim_results_tmp[1]==2 & n_case_FH < n_per_group) {
    sim_results <- cbind(sim_results, sim_results_tmp)
    n_case_FH <- n_case_FH+1
    #kept <- c(kept, 1) # for flagged info

  } 
  
  rm(sim_results_tmp)
  seed <- seed+1
}
#summary_indiv$kept <- kept # for flagged info


sim_results <- as.data.table(sim_results)
rownames(sim_results)[1:6]<-c("ordinal_score",
                              "disease", 
                              "num_family_affected",
                              "yri", 
                              "cconfounder", 
                              "seed")

#ordinal_score <- as.vector(t(sim_results[1, ]))
#tabulateordinal <- c(sum(ordinal_score == 0), sum(ordinal_score == 1), sum(ordinal_score == 2))
#text_to_save <-
#  paste0(
#    "CCC, ", paste(tabulateordinal, collapse = ", "), "\n",
#    "CCCpercent, ", paste0( formatC( tabulateordinal / ncol(sim_results), digits = 5, format = "f"), collapse = ","), "\n",
#    "variantstats_control, ", variant_statistics_freq_DT( sim_results[7:nrow(sim_results), which(ordinal_score == 0) , with = F], which_rare_causal), "\n",
#    "variantstats_caseNFH, ", variant_statistics_freq_DT( sim_results[7:nrow(sim_results), which(ordinal_score == 1) , with = F], which_rare_causal), "\n",
#    "variantstats_caseFH, ", variant_statistics_freq_DT( sim_results[7:nrow(sim_results), which(ordinal_score == 2) , with = F], which_rare_causal), "\n",
#    collapse = ", ")

## save simulation variant statistics
#write(text_to_save, file = paste0("variantstats_", outfile))

## save num_family_affected results
#t_results <- data.frame(t(sim_results[1:6,]))
#colnames(t_results) <- c("ordinal_score",
#                              "disease", 
#                              "num_family_affected",
#                              "yri", 
#                              "cconfounder", 
#                              "seed")
#num_rel_aff <- t_results$num_family_affected[t_results$ordinal_score==2]
#write.table(num_rel_aff,file=paste0("num_rel_affected_", outfile), row.names = FALSE, sep="\t", col.names = FALSE)

##------------------------------------------------------------------------------
## evaluate power of burden, skat, and gamut tests
##------------------------------------------------------------------------------
pvalues <- method_comparisons(sim_results=sim_results, residual_format=residual_format)
runtime <- unname((proc.time()-starttime)[3])


##------------------------------------------------------------------------------
## output results
##------------------------------------------------------------------------------
## fixit: write file row by row inside the for
df <- data.frame(pvalue.burden=pvalues[1], pvalue.skat=pvalues[2], pvalue.gamut=pvalues[3],
                 runtime=runtime, recomb.seed=recomb.seed, coal.seed=coal.seed, r.seed=r.seed)
write.table(df, outfile, row.names=FALSE, col.names=TRUE, sep="\t")


##------------------------------------------------------------------------------
## remove cosi .pos and .hap files
##------------------------------------------------------------------------------
file.remove(paste0("out",irep,".pos-1"),
            paste0("out",irep,".hap-1"),
            paste0("model",irep,".test"),
            paste0("recParams",irep),
            paste0("params",irep))
