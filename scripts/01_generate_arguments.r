## This script is run before the array of jobs 
## and saves all the arguments in logfile.txt and arguments.Rda

source("00_child_functions.r")

# simulation parameters
# =====================

# cosi parameters ** ONLY ONE ALLOWED **
seqlength <- 10000
nhap_eur <- 10000
nhap_afr <- 0

# simulation features
n_per_group <- c(750,1000)        # minimum number of people / case-control group

# family parameters
Nsibsize <- 0 #c(0, 1, 2, 3, 4)       # size of each sibship
minaff <- 1 # 2                      # minimum number of affected family members to be classified as "family history"
include_parents <- 1 # c(0, 1)         # are parents examined as part of "family history"?
# note: Nsibsize_eff <- Nsibsize + include_parents*2
randomize_nsibs <- 0               # c(0, 1)?

# rare variant effects
rare_maf <- 0.03                # rare-variant threshold
p_causal <- 0.15     # percent of rare variants that are causal
p_adverse <- 1 #c(0.5, 1)       # proportion of causal variants that increase risk
effectsize <- 1 #c(2, 4, 6)        # effect size of causal variant = log(4)/4
# note: beta_geno <- log(effectsize)/4
prev_target <- c(0.01,0.05) #c(0.05, 0.25)            # target prevalence
lambda <- c(2,4,8) #c(2,4,8)                          # conditional recurrence risk ratio

# covariate simulations
p_yri <- 0 # c(0, 0.25, 0.5)                 # population stratification
beta_yri <- 0 # 0.75 # 0                         # effect of YRI
cconfounder_family <- 0 #  c(0, 1)                 # continuous confounder family-aggregated?
beta_cconfounder <- 0.5 # 0.5                      # effect continuous confounder - change this to some % of variation

residual_format <- "PO" # "CR"

parameter_space <- expand.grid(seqlength, nhap_eur, nhap_afr, n_per_group,
                               Nsibsize, minaff, include_parents, randomize_nsibs,
                               rare_maf, p_causal, p_adverse, effectsize, prev_target, lambda,
                               p_yri, beta_yri, cconfounder_family, beta_cconfounder, residual_format)

names(parameter_space) <- c("seqlength", "nhap_eur", "nhap_afr", "n_per_group",
                           "Nsibsize", "minaff", "include_parents", "randomize_nsibs",
                           "rare_maf", "p_causal", "p_adverse", "effectsize", "prev_target", "lambda",
                           "p_yri", "beta_yri", "cconfounder_family", "beta_cconfounder", "residual_format")
parameter_space$residual_format <- as.character(parameter_space$residual_format)

# make sure to filter out....
# Nsibsize == 0 + include_parents == 0 
# also considering filter out high sibsize & high effect size (the cases w/o family history too rare)
# parameter_space <- parameter_space[ (parameter_space$Nsibsize != 0 & parameter_space$include_parents != 0), ] 

#rm(text,logfile)
#save.image(file="arguments.Rda")
save(parameter_space, file = "parameter_matrix.Rdata")
