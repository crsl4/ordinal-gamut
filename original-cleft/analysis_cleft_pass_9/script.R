dir <- "/Volumes/JFK_Passport/ordinal_family_history/analysis_cleft_pass_9/"
source(paste0(dir,"applied_functions.R"))

## study population files
load("/Volumes/JFK_Passport/ordinal_family_history/data_cleft/OFC_data_May2018.RData")

table(OFC$pcgroup)
table(OFC$AnalysisGroup_CaseControlOnly,OFC$pcgroup)
table(OFC$AnalysisGroup_CaseControlOnly,OFC$clefttype)
OFC$clefttype[OFC$clefttype==-9999] <- NA

## cleft lip (+ cleft lip&palette)

dat <- OFC[OFC$pcgroup %in% "CAUCASIAN",]
table(dat$AnalysisGroup_CaseControlOnly, dat$clefttype,useNA = "always")
dat$cleft_type_revised <- dat$clefttype
dat$cleft_type_revised[dat$cleft_type_revised %in% c(2,3)] <- 2.5
table(dat$AnalysisGroup_CaseControlOnly, dat$cleft_type_revised,useNA = "always")
dat$cleft_type_revised <- factor(dat$cleft_type_revised, levels=c(0,1,2.5,4),labels=c("unaffected",
                                                                                      "cleft but unknown type",
                                                                                      "cleft lip or lip+palate",
                                                                                      "cleft palate only"))
table(dat$AnalysisGroup_CaseControlOnly,dat$cleft_type_revised)
dat <- dat[dat$AnalysisGroup_CaseControlOnly%in%c("CASE","CONTROL"),]
table(dat$AnalysisGroup_CaseControlOnly,dat$cleft_type_revised)

dat_control <- dat[dat$AnalysisGroup_CaseControlOnly %in% "CONTROL",]
dat <- dat[dat$cleft_type_revised %in% "cleft lip or lip+palate",]
dat <- rbind(dat,dat_control)
table(dat$AnalysisGroup_CaseControlOnly,dat$cleft_type_revised)
dat$CLEFT <- 1
dat$CLEFT[dat$cleft_type_revised %in% "cleft lip or lip+palate"] <- 2
table(dat$CLEFT, dat$Multiplex, useNA = "ifany")

length(unique(dat$fid))
length(unique(paste(dat$fid,dat$iid,serp="/")))

## write ID list file to filter study population
ids_to_keep <- dat[,c("fid","iid")]
write.table(ids_to_keep,file=paste0(dir,"ids_to_keep.txt"),sep="\t",row.names = F,quote=F)

## write pheno file 
phenos_of_kept <- dat[,c("fid","iid","CLEFT")] # c("Multiplex","CAUC_EV1","CAUC_EV2","CAUC_EV3","CAUC_EV4","CAUC_EV5")
colnames(phenos_of_kept) <- c("FID","IID","CLEFT")
#phenos_of_kept$CLEFT[phenos_of_kept$CLEFT %in% "CASE"] <- 2
#phenos_of_kept$CLEFT[phenos_of_kept$CLEFT %in% "CONTROL"] <- 1
write.table(phenos_of_kept, file=paste0(dir,"phenos_of_kept.txt"),sep="\t",row.names = F,quote=F)

## run terminal_script

## write .cov file
fam <- read.delim(paste0(dir,"first_call/filtered.fam"), header=FALSE, stringsAsFactors=FALSE,sep=" ")
colnames(fam) <- c("FID","IID","PatIID","MatIID","Sex","Phenotype")

dat_ids <- paste(dat$fid,dat$iid,sep=":")
fam_ids <- paste(fam$FID,fam$IID,sep=":")

cov_text <- sapply(fam_ids, write_covs)
cov_text <- paste0(cov_text,collapse="")
cov_text <- paste0("FID IID FH PC1 PC2 PC3 PC4 PC5\n",cov_text)
write(cov_text,paste0(dir,"first_call/filtered.cov"))

## run bystro annotation

## format bystro annotation results
fn <- list.files(paste0(dir,"bystro"),pattern=".annotation.tsv")[1]
dat <- read.delim(paste0(dir,"bystro/",fn),stringsAsFactors = F)
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
write.table(gene_sets,paste0(dir,"bystro/bystro_gene_sets.txt"),sep="\t",row.names = F)

## create variant set ID file
gdat <- read.delim(paste0(dir,"bystro/bystro_gene_sets.txt"),sep="\t",stringsAsFactors = F)
g_uniq <- unique(gdat$gene)
logtxt <- paste0(logtxt, paste0('number of unique gene matches = ',length(g_uniq),"\n")) ## logtxt entry
bim <- read.delim(paste0(dir,"first_call/filtered.bim"), header=FALSE, stringsAsFactors=FALSE)
logtxt <- paste0(logtxt, paste0("number of plink-filtered variants used for annotation = ",nrow(bim),"\n")) ## logtxt entry
colnames(bim) <- c("chrom","vid","pos_morg","coord","allele1","allele2")
bim$marker <- paste(bim$chrom,bim$coord,sep=":")

min_rv_per_gene=4
set_list <- sapply(g_uniq, write_gene_variants, min_rv_per_gene=min_rv_per_gene)
good_genes <- names(set_list)[!set_list==""]
length(good_genes)
logtxt <- paste0(logtxt, paste0("min RV per gene filter = ",min_rv_per_gene),"\n") ## logtxt entry
set_list <- paste0(set_list, collapse="")
write(set_list,paste0(dir,"good_genes.set"))

## create by_gene folders and gene-level plink files
call <- sapply(good_genes, by_gene_plink, dir)
call <- paste0(call, collapse="")
write(call, paste0(dir,"set_call_script"))

## run set_call_script

### analyze sets
pheno_mat <- "linear kernel"
cov <- read.delim(paste0(dir,"first_call/filtered.cov"),sep=" ",stringsAsFactors = F)
rev_good_genes <- list.dirs(paste0(dir,"by_gene"), full.names = F)[-1] ## why missing some folders?
logtxt <- paste0(logtxt, paste0("number of tests (num genes meeting min RV threshold) = ",length(rev_good_genes),"\n")) ## logtxt entry
logtxt <- paste0(logtxt, paste0("phenotypic similarity = ",pheno_mat),"\n") ## logtxt entry
write(logtxt,file=paste0(dir,"logfile.txt"))

#which(rev_good_genes%in%"C1orf112")
rev_good_genes_2 <- read.delim(paste0(dir,"by_gene/gene_test_results.txt"),stringsAsFactors = F, sep=" ",header = F)
rev_good_genes_2[nrow(rev_good_genes_2),]
rev_good_genes[2366]
rev_good_genes <- rev_good_genes[-1:-2366]
head(rev_good_genes)
results <- sapply(rev_good_genes, do_gene_test, cov=cov, dir=dir, pheno_mat=pheno_mat)

## read in results
res <- read.delim(paste0(dir,"by_gene/gene_test_results.txt"),stringsAsFactors = F, sep=" ",header = F)
colnames(res) <- c("gene","p.burden","p.skat","p.gamut")
gdat <- read.delim(paste0(dir,"bystro/bystro_gene_sets.txt"),sep="\t",stringsAsFactors = F)

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
save(res_with_loc,file=paste0(dir,"by_gene/res_with_loc.Rda"))

 