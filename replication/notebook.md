# Replication dataset

## Data

Data downloaded from Elizabeth Leslie's email: `20181212_data_us.tar.gz`. This folder contains:
```shell
$ cd Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/
$ ls
20181212_data_us.tar.gz data_us
$ cd data_us/
$ ls
data_recode_8:101600001-106200000.bed data_recode_8:101600001-106200000.fam
data_recode_8:101600001-106200000.bim data_recode_8:101600001-106200000.log
```

Info from emails:
- Elizabeth: We would need a raw phenotype definition (probably cleft type), family history yes/no if you have it, and raw genetic data. We would only need the region for 8q22.3 (hg19 chr8:101,600,001-106,200,000) , so that should be easier than the whole dataset.
- In the attachment, you will find files prepared by Kacper. He told me that he excluded SNPs with SNP call rate <0.95 and samples with sample call rate < 0.95.


## Analysis

We will follow Taylor's script used for the original analysis.
The folder `analysis_cleft_pass_9` has the scripts from Taylor's analysis (see `readme.md`), and `results_cleft` has the Rmd script for results.

```shell
cd Dropbox/Documents/gwas/projects/ordinal-skat/replication/
scp csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/ordinal_family_history_files_for_claudia/analysis_cleft_pass_9/* analysis_cleft_pass_9
scp csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/ordinal_family_history_files_for_claudia/results_cleft/* results_cleft
```

In the original analysis that Taylor ran, she used 4 input files:
```shell
hgcc:node00:[data_cleft] % ls
CIDR_Marazita_OFC_TOP_subject_level_filtered.bed
CIDR_Marazita_OFC_TOP_subject_level_filtered.bim
CIDR_Marazita_OFC_TOP_subject_level_filtered.fam
OFC_data_May2018.RData
```
According to the R scripts, the RData file has the phenotype, which was later added to the bedfiles with plink.

For the replication analysis, we received only 3 files:
```shell
$ ls
data_recode_8:101600001-106200000.bed 
data_recode_8:101600001-106200000.fam
data_recode_8:101600001-106200000.bim 
```

I am assuming that these files already have the phenotype that we want to study.
Based on email from Elizabeth, it seems that this is the correct phenotype. We might want to check the values:
```shell
$ head *.fam
1 CLP1 0 0 1 2
1 CLP10 0 0 1 2
1 CLP100 0 0 1 2
1 CLP101 0 0 1 2
1 CLP102 0 0 1 2
1 CLP103 0 0 1 2
1 CLP104 0 0 2 2
1 CLP105 0 0 1 2
1 CLP106 0 0 2 2
1 CLP107 0 0 2 2

$ wc -l data_recode_8\:101600001-106200000.fam
     845 data_recode_8:101600001-106200000.fam
```

We can read the fam file in R and study the 6th column. 

```r
setwd("Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us")
dat = read.table("data_recode_8:101600001-106200000.fam", header=FALSE, sep=' ')

> summary(as.factor(dat$V6))
  1   2 
570 275
```

We only have two categories (1=controls, 2=cases), instead of 3 as in the original analysis.
Also, we do not have family information, so we will need to run a usual SKAT.

After talking to Mike, we want to check the common variants first.

1. Standard GWAS common variants

First, we only consider common variants:
```shell
plink --bfile data_recode_8:101600001-106200000 --maf 0.05 --make-bed --out common-variants

8192 MB RAM detected; reserving 4096 MB for main workspace.
1204 variants loaded from .bim file.
845 people (438 males, 407 females) loaded from .fam.
845 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 845 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997705.
169 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
1035 variants and 845 people pass filters and QC.
Among remaining phenotypes, 275 are cases and 570 are controls.
--make-bed to common-variants.bed + common-variants.bim + common-variants.fam
... done.

$ wc -l common-variants.bim
    1035 common-variants.bim
```

But, it turns out that there are some repeated snps:
```shell
$ grep "exm-" common-variants.bim
8	exm-rs3108919	0	101841546	G	A
8	exm-rs1264202	0	102227923	A	C
8	exm-rs2458413	0	105359432	G	A
8	exm-rs16872085	0	105957540	G	A
8	exm-rs977396	0	105992696	G	A
8	exm-rs9297357	0	106142333	G	A
```
Because for example exm-rs3108919 and rs3108919 both exist, and they are both complete.

So, we will remove these from the analysis:
```shell
grep "exm-" common-variants.bim | awk '{ print $2 }' > exclude-exm.txt
```

And we exclude those snps:
```shell
plink --bfile common-variants --exclude exclude-exm.txt --make-bed --out common-variants2
```

Now, we run association study with plink on these common variants:
```shell
plink --bfile common-variants2 --assoc
```

Now, we look at the pvalues to see if there is any inflation (copied `qqplot-MF.r`).
```r
source("../../../scripts/qqplot-MF.r")
dat = read.table("plink.assoc", header=TRUE)
pdf("qqplot-common-variants.pdf", width=5, height=5)
qqunif.plot(dat$P)
dev.off()
```
Based on QQplot, there is no evidence of inflation, so we will not adjust by PCs here.

Now, we want to see if there is anything significant with Bonferroni correction:
```r
ntests = length(dat$P)
alpha = 0.05/ntests
str(dat[dat$P < alpha,]) ##nothing significant

str(dat[dat$P < 0.0005,])

> dat[dat$P < 0.0005,]
    CHR      SNP        BP A1   F_A    F_U A2 CHISQ         P    OR
491   8 rs965660 103597492  G 0.173 0.1004  A 14.52 0.0001389 1.874
```

Now, we want to include sex as covariate:
```shell
plink --bfile common-variants2 --logistic --sex --out plink-sex
```

Now, we look at the pvalues to see if there is any inflation (copied `qqplot-MF.r`).
```r
source("../../../scripts/qqplot-MF.r")
dat = read.table("plink-sex.assoc.logistic", header=TRUE)
dat2 = subset(dat,TEST=="ADD") ##removing pvalues for SEX
pdf("qqplot-common-variants-sex.pdf", width=5, height=5)
qqunif.plot(dat2$P)
dev.off()
```

Based on QQplot, there is no evidence of inflation, so we will not adjust by PCs here.

Now, we want to see if there is anything significant with Bonferroni correction:
```r
ntests = length(dat2$P)
alpha = 0.05/ntests
str(dat2[dat2$P < alpha,]) ##nothing significant

str(dat2[dat2$P < 0.0005,])

> dat2[dat2$P < 0.0005,]
    CHR      SNP        BP A1 TEST NMISS    OR STAT         P
981   8 rs965660 103597492  G  ADD   704 2.014 3.89 0.0001001
```

According to Mike, we do not want to do Bonferroni correction because the tests are not independent. So, we will use `p_ACT` R script (from Karen Conneely) to adjust for multiple testing on correlated test.

### p_ACT for multiple testing

1. Downloaded `p_ACT` from [here](http://csg.sph.umich.edu/boehnke/p_act/p_ACT_1.2/p_act.php).

2. Created `alpha.txt` file with 0.05 as significance level in `scripts/p_ACT_1.2`.

3. We need to read `plink-sex.assoc.logistic` and create the `pvalues.txt` with 3 columns: trait, SNP, pvalue. We do this in R. The file needs to be space-delimited:
```r
setwd("Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us")
dat = read.table("plink-sex.assoc.logistic", header=TRUE)
dat2 = subset(dat,TEST=="ADD") ##removing pvalues for SEX
dat3 = data.frame(y="y", snp=dat2$SNP, pval=dat2$P)

> head(dat3)
  y        snp   pval
1 y  exm713146 0.5788
2 y  rs2096548 0.7014
3 y  rs2469666 0.5880
4 y   rs752036 0.6094
5 y rs10111792 0.6808
6 y  rs7017809 0.1752

## Remove 35 with higher pvalue
q = unname(quantile(dat3$pval,1000/length(dat3$pval)))
dat4 = dat3[dat3$pval <= q,]

write.table(dat4, file="pvalues.txt", sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)
```
Note that there is a hard limit for `mvtnorm` of 1000 tests. So, I need to get rid of 35.

4. Now we need the `genotype.txt` file with the genotypes for all 1035 snps. Currently, we have the files `common-variants2.bed/bim/fam`, we need to use plink to write as "raw data".
```shell
plink --bfile common-variants2 --recodeA --out common-variants-raw
```
This creates the `common-variants-raw.raw` which has the genotypes 0/1/2 with column names given by the specific SNPs. This is almost what we need, except that the column names have the allele added. For example, the first SNP is `exm713146_G` (column name in `common-variants-raw.raw`), but it is named `exm713146` in the `pvalues.txt` file.

This raw file also has extra columns: `FID IID PAT MAT SEX PHENOTYPE` not needed in `genotype.txt`.

So, we read in R to modify:
```r
dat = read.table("common-variants-raw.raw", header=TRUE)
dat$FID = NULL
dat$IID = NULL
dat$PAT = NULL
dat$MAT = NULL
dat$SEX = NULL
dat$PHENOTYPE = NULL

for(i in 1:length(colnames(dat))){
    s = strsplit(colnames(dat)[i], "_")
    colnames(dat)[i] = s[[1]][1]    
}

## Now we remove the 29 with worse pvalue
pval = read.table("pvalues.txt", header=FALSE)
bad.snps = setdiff(colnames(dat), as.character(pval$V2))
dat2 = dat[ , !(names(dat) %in% bad.snps)]

## they are kept in same order:
> sum(colnames(dat2) == as.character(pval$V2))
[1] 1000

write.table(dat2, file="genotype.txt", row.names=FALSE, sep=' ', quote=FALSE)
```

5. Now we have in `data/data_us` the files: `alpha.txt, genotype.txt, pvalues.txt`, and we will run `p_ACT_seq.R`. We tried to run in the terminal, but we got an error:
```shell
$ Rscript ../../scripts/p_ACT_1.2/p_ACT_seq.R 
Warning message:
package ‘mvtnorm’ was built under R version 3.4.4 
Read 1 item
## it takes seconds!

$ ls *out
p_ACT_seq.out

$ head p_ACT_seq.out 
p_ACT_seq.R, version 1.2 
************************************************
Sequential p_ACT for smallest p-values 
 
Trait	       SNP	        p	   p_ACT
    y	  rs965660	0.0001001	0.064419
    y	rs10094493	0.0012890	      NA
    y	 rs2511669	0.0027420	      NA
    y	rs10955341	0.0048480	      NA
    y	rs11785817	0.0079980	      NA
```

We can also try the non-sequential version for fun:
```shell
$ Rscript ../../scripts/p_ACT_1.2/p_ACT.R 
Warning message:
package ‘mvtnorm’ was built under R version 3.4.4

$ ls *out
p_ACT.out     p_ACT_seq.out

$ head p_ACT.out 
p_ACT.R, version 1.2 
************************************************
1000 tests performed. 
Trait and SNP with smallest p-value:  y rs965660 
Smallest p-value:  0.0001001 
p_ACT: 1 
```

Thus, no significant hit in common variants!


## Following Taylor's script

1. Now, we want to follow `script.R` for the tests of rare variants. First, we want to extract only rare variants based on `terminal_script`:
```shell
cd Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us
plink --bfile data_recode_8:101600001-106200000 --maf 0.001 --max-maf 0.05 --missing --freq --make-bed --out filtered

65536 MB RAM detected; reserving 32768 MB for main workspace.
1204 variants loaded from .bim file.
845 people (438 males, 407 females) loaded from .fam.
845 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 845 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997705.
--freq: Allele frequencies (founders only) written to filtered.frq .
--missing: Sample missing data report written to filtered.imiss, and
variant-based missing data report written to filtered.lmiss.
1068 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
136 variants and 845 people pass filters and QC.
Among remaining phenotypes, 275 are cases and 570 are controls.
--make-bed to filtered.bed + filtered.bim + filtered.fam ... done.

$ wc -l filtered.bim 
     136 filtered.bim
```
We manually checked if there are any of the "exm-" markers, but no.

2. Now, we want to check if all these snps are already in the bystro output file.

```r
bim = read.table("filtered.bim", header=FALSE)
goodgenes = read.table("../../analysis_cleft_pass_9/good_genes.set", header=FALSE)

snp=rep(NA, length(bim$V2))
present= rep(NA,length(bim$V2))
for(i in 1:length(bim$V2)){
    s = sum(as.character(goodgenes$V1) == as.character(bim$V2[i]))
    snp[i] = as.character(bim$V2[i])
    present[i] = s
}

df = data.frame(snp=snp, present=present)

> sum(df$present)
[1] 13

write.table(df,file="shared-snps.csv", row.names=FALSE, quote=FALSE)
```
We share only 13 markers with the original data.


3. Create `vcf_filtered` files in plink

Following the `script.R`, we do not need to extract given individuals or SNPs, so we only run the second command of `terminal_script`:
```shell
cd Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us
plink --bfile filtered --recode vcf --out vcf_filtered 

PLINK v1.90b4.4 64-bit (21 May 2017)           www.cog-genomics.org/plink/1.9/
(C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to vcf_filtered.log.
Options in effect:
  --bfile filtered
  --out vcf_filtered
  --recode vcf

65536 MB RAM detected; reserving 32768 MB for main workspace.
136 variants loaded from .bim file.
845 people (438 males, 407 females) loaded from .fam.
845 phenotype values loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 845 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.99826.
136 variants and 845 people pass filters and QC.
Among remaining phenotypes, 275 are cases and 570 are controls.
--recode vcf to vcf_filtered.vcf ... done.
```

4. Running bystro. Based on Taylor's email (`readme.md`), we ran all the steps and downloaded a file `vcf_filtered.tar`, which we moved to the data folder. After untarring the files, I have `vcf_filtered.annotation.tsv`.

5. Format bystro annotation results (based on `script.R`):
```r
setwd("Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us/vcf_filtered")
fn <- "vcf_filtered.annotation.tsv"
dat <- read.delim(fn,stringsAsFactors = F)
logtxt <- paste0("number of bystro annotations = ", nrow(dat),"\n") ## logtxt start
set <- data.frame(matrix(nrow=0,ncol=5))
colnames(set) <- c("chrom","pos","gnomad.genomes.id","type","gene")

for(i in 1:nrow(dat)){
  tmp <- data.frame(matrix(nrow=0,ncol=5))
  colnames(tmp) <- c("chrom","pos","gnomad.genomes.id","type","gene")
  
  tmp[i,c("chrom","pos","gnomad.genomes.id")] <- dat[i,c("chrom","pos","gnomad.genomes.id")]
  genes <- unique(unlist(strsplit(dat[i,"refSeq.name2"],";")))
  nearest_genes <- unique(unlist(strsplit(dat[i,"refSeq.nearest.name2"],";")))
  tmp <- tmp[rep(i,length(genes)+length(nearest_genes)),]
  tmp$type <- c(rep("refSeq.name2",length(genes)),rep("refSeq.nearest.name2",length(nearest_genes)))
  tmp$gene <- c(genes, nearest_genes)
  set <- rbind(set,tmp)
  rm(tmp); rm(genes); rm(nearest_genes)
  print(paste0("row ",i," completed"))
}

gene_sets <- set[set$type%in%"refSeq.name2",]
unique_genes <- unique(gene_sets$gene)
length(unique_genes)
gene_sets$marker <-  paste(substr(gene_sets$chrom,4,nchar(gene_sets$chrom)),gene_sets$pos,sep=":")
write.table(gene_sets,file="bystro_gene_sets.txt",sep="\t",row.names = F)

## create variant set ID file
gdat <- read.delim("bystro_gene_sets.txt",sep="\t",stringsAsFactors = F)
g_uniq <- unique(gdat$gene)
logtxt <- paste0(logtxt, paste0('number of unique gene matches = ',length(g_uniq),"\n")) ## logtxt entry
bim <- read.delim("../filtered.bim", header=FALSE, stringsAsFactors=FALSE)
logtxt <- paste0(logtxt, paste0("number of plink-filtered variants used for annotation = ",nrow(bim),"\n")) ## logtxt entry
colnames(bim) <- c("chrom","vid","pos_morg","coord","allele1","allele2")
bim$marker <- paste(bim$chrom,bim$coord,sep=":")

min_rv_per_gene=4
source("../../../analysis_cleft_pass_9/applied_functions.R")
set_list <- sapply(g_uniq, write_gene_variants, min_rv_per_gene=min_rv_per_gene)
good_genes <- names(set_list)[!set_list==""]
length(good_genes) ##9
logtxt <- paste0(logtxt, paste0("min RV per gene filter = ",min_rv_per_gene),"\n") ## logtxt entry
set_list <- paste0(set_list, collapse="")
write(set_list,"good_genes.set")
write(logtxt,file="logfile-bystro.txt")

## create by_gene folders and gene-level plink files
dir=""
call <- sapply(good_genes, by_gene_plink, dir)
call <- paste0(call, collapse="")
write(call, "set_call_script")
```

6. **not reproducible step: manually done** We will copy `set_call_script` and `good_genes.set` outside of `vcf_filtered` bystro folder, because we will modify manually. First, we will get rid of gene denoted ! (which represents missing info in bystro).
We also need to get rid of `../../` before `plink`, and remove the `../../first_call/` and change to `../../`. We could have modified the `by_gene_plink` function instead, but did not do this.

7. We create a `by_gene` folder inside `data/data_us`. We copy and paste lines in `set_call_script`.

8. Now, we need to create the `filtered.cov` file with the family history (NA here) and the PCs. However, I do not understand where this is created, so I will copy the files from hgcc:
```shell
cd Dropbox/Documents/gwas/projects/ordinal-skat/replication/analysis_cleft_pass_9/first_call
scp csolislemus@hgcc.genetics.emory.edu:/home/csolislemus/ordinal_family_history_files_for_claudia/analysis_cleft_pass_9/first_call/filtered* .
```
It turns out that the PCs were already given to Taylor in `OFC_data_May2018.RData`, and we do not have them. I will ignore the PCs for now (based on common variants qqplot).

9. Now, we run the gene tests (based on `script.R`). We need to modify the do_gene_test function to only run burden and skat, and also, so that it does not use the `cov` file (with FH and PCs, that we do not have for the replication dataset)
```r
do_gene_test <- function(gene) {
  
  file <- paste0("by_gene/",gene,"/setdat.raw")
  raw <- mean_impute(file)
  
  
  ##FH <- cov$FH ## assumes FH is 0/1 variable and first covariates listed in covar file
  ##COVAR <- as.matrix(cov[,-1:-3])
  PHENO <- (raw$PHENOTYPE)-1
  GENO <- as.matrix(raw[,-1:-6])
  
  ## 0 = control, 1 = case w/o family history, 2 = case w/ history
  ##ordinal_score <- PHENO + PHENO*(FH) # family history
  
  # remove monomorphic SNPs
  GENO <- t(GENO) # rows are SNPs and columns are subjects
  GENO <- GENO[rowSums(GENO) != 0, ]
  
  # weighting (beta weights via original SKAT)
  MAF_estimate = rowMeans(GENO)/2
  beta_weight = dbeta(MAF_estimate, 1, 25) / dbeta(0, 1, 25) 
  
  # burden
  geno_burden <- colSums(beta_weight %*% GENO) # weighted sum
  model <- glm(PHENO ~ geno_burden, family = binomial) 
  burden <- summary(model)$coefficients[2, 1]
  
  # skat
  skat_obj <- SKAT_Null_Model(PHENO ~ 1, out_type = "D") 
  skat <- SKAT(Z = t(GENO), obj = skat_obj, is_check_genotype = F, weights = beta_weight)$p.value
  
  p <- c(burden, skat)
  p <- paste(gene, paste0(p,collapse=" "))
  write(p, "by_gene/gene_test_results.txt", append = T)
  print(p)
  
  return(p)
}
```

And we run the gene tests:
```r
setwd("Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us")
source("../../analysis_cleft_pass_9/applied_functions.R")
rev_good_genes <- list.dirs("by_gene", full.names = F)[-1]
results <- sapply(rev_good_genes, do_gene_test)
```

10. Now, we summarize the results
```r
## read in results
setwd("Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us")
res <- read.delim("by_gene/gene_test_results.txt",stringsAsFactors = F, sep=" ",header = F)
colnames(res) <- c("gene","p.burden","p.skat")
gdat <- read.delim("vcf_filtered/bystro_gene_sets.txt",sep="\t",stringsAsFactors = F)

loc <- sapply(res$gene,function(x) {
  tmp <- unique(gdat$chrom[which(gdat$gene%in%x)])[1] 
  chrom <- as.numeric(substr(tmp,4,nchar(tmp)))
  tmp <-  unique(gdat$pos[which(gdat$gene%in%x)])[1] 
  return(chrom)
}
)

pos <- sapply(res$gene,function(x) {
  pos <-  unique(gdat$pos[which(gdat$gene%in%x)])[1] 
  return(pos)
}
)

res_with_loc <- cbind(res,loc,pos)
save(res_with_loc,file="by_gene/res_with_loc.Rda")
```

11. Summary plots (from `cleft_results_pass9.Rmd`)

```r
setwd("Dropbox/Documents/gwas/projects/ordinal-skat/replication/data/data_us")
##source("../../analysis_cleft_pass_9/applied_functions.R")
load("by_gene/res_with_loc.Rda")
colnames(res_with_loc) <- c("gene","p.burden","p.skat","chr","bp")
paste0("number of tests (ie. number of genes with at least 4 rare variants) = ",nrow(res_with_loc))

## significance thresholds
sugg = -log10(1e-04) ##=4
gen = -log10(0.05/nrow(res_with_loc)) ##= 2.20412 

## Plots
library(ggplot2)
library(qqman)

pdf("replication-cleft.pdf")
manhattan(res_with_loc, chr = "chr", bp = "bp", p = "p.skat",snp="gene",
  col = c("blue4", "orange3"), chrlabs = NULL, annotatePval = 1e-3,
  suggestiveline = sugg, genomewideline = gen, main="SKAT",ylim=c(0,6), xlim=c(9e7, 1.2e8))

qq(res_with_loc$p.skat, main = "Q-Q plot of SKAT p-values")

manhattan(res_with_loc, chr = "chr", bp = "bp", p = "p.burden",snp="gene",
  col = c("blue4", "orange3"), chrlabs = NULL, annotatePval = 1e-3,
  suggestiveline = sugg, genomewideline = gen, main="burden",ylim=c(0,6), xlim=c(9e7, 1.2e8))

qq(res_with_loc$p.burden, main = "Q-Q plot of burden test p-values")
dev.off()
```

Nothing significant.