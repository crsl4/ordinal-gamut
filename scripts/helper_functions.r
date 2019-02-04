
##------------------------------------------------------------------------------
## evaluate type 1 error
##------------------------------------------------------------------------------
eval_type1_error <- function (pvals, threshold) {
  apply(pvals, 2, function(x) { sum(x <= threshold) / length(x) })
}


##------------------------------------------------------------------------------
## split up runs by a parameter
##------------------------------------------------------------------------------
split_runs <- function(sim_folder, split_by) {
  list_of_folders <- list.files(sim_folder, pattern = "^run_")
  load(paste0(sim_folder,"parameter_matrix.Rdata"))
  levs <- unique(parameter_space[,split_by])
  levlist <- list()
  for(i in 1:length(levs)){
    levlist[[i]] <- which(parameter_space[,split_by]%in%levs[i])
    names(levlist)[[i]] <- levs[i]
  }
  return(list(sim_folder=sim_folder,levlist=levlist,split_by=split_by))
}


##------------------------------------------------------------------------------
## read in one file of run results
##------------------------------------------------------------------------------
read_one_sim <- function(file) {
  dat <- read.delim(file, stringsAsFactors=FALSE)
  tmp <- unlist(strsplit(file,"ordinal_simulation_"))[2]
  dat$irep <- as.numeric(substr(tmp,1,nchar(tmp)-4))
  dat$filepath <- file
  return(dat=dat)
}


##------------------------------------------------------------------------------
## aggregate and save simulation results
## sim_folder needs trailing slash /
## output_dir does not need trailing slash
##------------------------------------------------------------------------------
#process_sim_files(sim_folder="/Volumes/JFK_Passport/ordinal_family_history/code/ordinal_flagged_201805241446_post/",
#                  output_dir="/Volumes/JFK_Passport/ordinal_family_history/code/ordinal_flagged_201805241446_post/processed_data")
process_sim_files <- function(sim_folder, filter_by=NA, output_dir) {
  library(data.table)
  library(reshape2)
  
  list_of_folders <- list.files(sim_folder, pattern = "^run_")
  load(paste0(sim_folder,"parameter_matrix.Rdata"))
  
  ## check
  length(list_of_folders)==nrow(parameter_space)
  
  if(!is.na(filter_by)){
    parameter_space <- parameter_space[filter_by,]
    list_of_folders <- list_of_folders[which(list_of_folders%in%paste0("run_",filter_by))]
  }
  
  list_of_folders <- list_of_folders[order(as.numeric(substr(list_of_folders,5,nchar(list_of_folders))))]
  
  res <- lapply(list_of_folders, function(x) {
    list_of_files <- list.files(paste0(sim_folder,x), pattern = "^ordinal_simulation_", full.names = TRUE)
    l <- lapply(list_of_files, read_one_sim)
    df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T),stringsAsFactors=FALSE)
    colnames(df) <- names(l[[1]])
    return(df)
  })
  
  names(res) <- list_of_folders
  for(i in 1:length(res)){
    run_table <- res[[i]]
    param_table <- parameter_space[i,]
    save(run_table, param_table, file=paste0(output_dir,"/results_",names(res)[i],".RData"))
  }
  
}

##------------------------------------------------------------------------------
## calculate power
##------------------------------------------------------------------------------
get_power <- function(file, alpha) {
  load(file)
  grab <- grep("pvalue",colnames(run_table))
  power <- apply(run_table[,grab],2,function(x) sum(as.numeric(x) < alpha)/length(x))
  tmp <- cbind(power,nrow(run_table),param_table)
  colnames(tmp)[2] <- "nsims"
  tmp$method <- unlist(strsplit(rownames(tmp),"pvalue."))[!unlist(strsplit(rownames(tmp),"pvalue."))==""]
  tmp <- tmp[,c(ncol(tmp),1:(ncol(tmp)-1))]
  row.names(tmp)<-NULL
  rm(run_table)
  rm(param_table)
  return(power_table=tmp)
}


##------------------------------------------------------------------------------
## plot power
##------------------------------------------------------------------------------
plot_power <- function(input_dir, alpha, xvar="n_per_group", xlabel="Samples/Group", topvar="Nsibsize", topprefix="N Siblings = ",
                       rightvar="effectsize", rightprefix="RV Effect Size = ", linevar="method", linelabel="Method", plot_ymax=1, filter_by=NA) {
  
  if(length(input_dir)>1) {
    list_of_files <- unlist(lapply(input_dir, list.files, full.names=TRUE))
  } else {
    list_of_files <- list.files(input_dir,full.names = TRUE)
  }
  
  l <- lapply(list_of_files, get_power, alpha=alpha)
  
  dat <- {}
  for(i in 1:length(l)) {
    dat <- rbind(dat,l[[i]])
  }
  
  if(!is.na(filter_by)) {
    for(i in 1:length(filter_by)){
      varid <- names(filter_by[i])
      vals <- filter_by[[i]]
      dat <- dat[dat[,varid]%in%vals,]
    }
  }
  
  constant_params <- dat[1,apply(dat, 2, function(col) { length(unique(col))==1 })]
  row.names(constant_params)<-"value"
  
  dat$method <- factor(dat$method,levels=c("burden","skat","gamut"),labels=c("Burden Test","SKAT","Ordinal GAMuT"))
  if(!linevar%in%"method") {
    dat[,linevar] <- factor(dat[,linevar])
  }
  dat[,topvar] <- factor(dat[,topvar])
  levels(dat[,topvar]) <- paste0(topprefix,levels(dat[,topvar]))
  dat[,rightvar] <- factor(dat[,rightvar])
  levels(dat[,rightvar]) <- paste0(rightprefix,levels(dat[,rightvar]))
  
  library(ggplot2)
    fig <- ggplot(data = dat, aes_string(x = xvar, y = "power", color = linevar)) + 
    geom_point(aes_string(shape = linevar)) + geom_line(aes_string(group = linevar), alpha = 0.6)  + 
    scale_x_continuous(labels = unique(dat[,xvar]), breaks = unique(dat[,xvar])) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(limits = c(0, plot_ymax), expand = c(0,0)) + 
    facet_grid(paste0(rightvar,"~",topvar)) + guides(color=guide_legend(title=linelabel),shape="none") +
    xlab(xlabel) + ylab(paste0("Empirical Power (α = ",alpha,")"))
    
    return(list(fig=fig, dat=dat, constant_params=constant_params))
}


##------------------------------------------------------------------------------
## read in one file of variant stat results
##------------------------------------------------------------------------------
read_one_variant <- function(file) {
  dat <- read.delim(file, stringsAsFactors = F, sep=",", header = FALSE)
  formed <- data.frame(matrix(nrow=3,ncol=7))
  colnames(formed)<-c("group","group_count","group_prop","avg_MAF","avg_n_variants","avg_n_causal","irep")
  formed$group <- c("control","caseNFH","caseFH")
  formed$group_count <- as.numeric(dat[1,-1])
  formed$group_prop <- as.numeric(dat[2,-1])
  formed$avg_MAF <- as.numeric(dat[-1:-2,2])
  formed$avg_n_variants <- as.numeric(dat[-1:-2,3])
  formed$avg_n_causal <- as.numeric(dat[-1:-2,4])
  
  tmp <- unlist(strsplit(file,"variantstats_ordinal_simulation_"))[2]
  formed$irep <- as.numeric(substr(tmp,1,nchar(tmp)-4))
  formed$filepath <- file
  return(dat=formed)
}


##------------------------------------------------------------------------------
## aggregate and save variant stats per group
## sim_folder does not need trailing slash
## output_dir does not need trailing slash
##------------------------------------------------------------------------------
process_variant_files <- function(sim_folder, output_dir, filter_by=NA) {
  library(data.table)
  library(reshape2)
  
  list_of_folders <- list.files(sim_folder, pattern = "^run_")
  load(paste0(sim_folder,"/parameter_matrix.Rdata"))
  
  ## check
  length(list_of_folders)==nrow(parameter_space)
  
  if(!is.na(filter_by)){
    parameter_space <- parameter_space[filter_by,]
    list_of_folders <- list_of_folders[which(list_of_folders%in%paste0("run_",filter_by))]
  }
  
  list_of_folders <- list_of_folders[order(as.numeric(substr(list_of_folders,5,nchar(list_of_folders))))]
  
  res <- lapply(list_of_folders, function(x) {
    list_of_files <- list.files(paste0(sim_folder,"/",x), pattern = "^variantstats_", full.names = TRUE)
    l <- lapply(list_of_files, read_one_variant)
    df <- {}
    for(i in 1:length(l)) {
      df <- rbind(df, l[[i]])
    }
    return(df)
  })
  
  names(res) <- list_of_folders
  for(i in 1:length(res)){
    run_table <- res[[i]]
    param_table <- parameter_space[i,]
    save(run_table, param_table, file=paste0(output_dir,"/results_",names(res)[i],".RData"))
  }
  
}


##------------------------------------------------------------------------------
## helper var stat 
##------------------------------------------------------------------------------
get_varstat <- function(file) {
  load(file)
  run_table$group <- factor(run_table$group, levels=c("control","caseNFH","caseFH"),labels=c("Control","Case No FH","Case FH"))
  dat <- data.frame(group=levels(run_table$group),
                    stat=rep("avg_MAF",3),
                    mean=tapply(run_table$avg_MAF,run_table$group,mean),
                    sd=tapply(run_table$avg_MAF,run_table$group,sd))
  dat <- cbind(dat, rbind(param_table,param_table,param_table))
  row.names(dat)<-NULL
  dat2 <- data.frame(group=levels(run_table$group),
                    stat=rep("avg_n_variants",3),
                    mean=tapply(run_table$avg_n_variants,run_table$group,mean),
                    sd=tapply(run_table$avg_n_variants,run_table$group,sd))
  dat2 <- cbind(dat2, rbind(param_table,param_table,param_table))
  row.names(dat2)<-NULL
  dat3 <- data.frame(group=levels(run_table$group),
                     stat=rep("avg_n_causal",3),
                     mean=tapply(run_table$avg_n_causal,run_table$group,mean),
                     sd=tapply(run_table$avg_n_causal,run_table$group,sd))
  dat3 <- cbind(dat3, rbind(param_table,param_table,param_table))
  row.names(dat3)<-NULL
  
  tmp <- rbind(dat,dat2,dat3)
  rm(run_table)
  rm(param_table)
  return(dat=tmp)
}


##------------------------------------------------------------------------------
## analyze variant stats per group
##------------------------------------------------------------------------------
#input_dir <- "/Volumes/JFK_Passport/ordinal_family_history/code/ordinal_201803081323_post/processed_variant_data"
analyze_varstats <- function(input_dir, filter_by=NA) {
  
  if(length(input_dir)>1) {
    list_of_files <- unlist(lapply(input_dir, list.files, full.names=TRUE))
  } else {
    list_of_files <- list.files(input_dir,full.names = TRUE)
  }
  
  l <- lapply(list_of_files, get_varstat)
  
  dat <- {}
  for(i in 1:length(l)) {
    dat <- rbind(dat,l[[i]])
  }
  
  if(!is.na(filter_by)) {
    for(i in 1:length(filter_by)){
      varid <- names(filter_by[i])
      vals <- filter_by[[i]]
      dat <- dat[dat[,varid]%in%vals,]
    }
  }
  
  library(ggplot2)
  
  dodge <- position_dodge(width = 0.9)
  limits <- aes(ymax = dat$mean + dat$sd,
                ymin = dat$mean - dat$sd)
  
  p <- ggplot(data = dat, aes(x = group, y = mean, fill = group)) + facet_grid(effectsize~Nsibsize)
  
  p + geom_bar(stat = "identity", position = dodge) +
    geom_errorbar(limits, position = dodge, width = 0.25) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.x=element_blank())
  
}


##------------------------------------------------------------------------------
## read in one file of flagged output 
##------------------------------------------------------------------------------
read_one_flagged <- function(file) {
  
  dat <- read.delim(file, stringsAsFactors = F, header = TRUE)
  
  dat$group <- paste0(dat$ordinal.score,"_kept=",dat$kept)
  dat$group <- factor(dat$group)
  
  if(length(grep("^pdis.sib",colnames(dat)))>1){
    dat$pdis.sib.avg <- apply(dat[,grep("^pdis.sib",colnames(dat))],1,mean)
    dat$rev.pdis.sib.avg <- apply(dat[,grep("^rev.pdis.sib",colnames(dat))],1,mean)
  } else {
    dat$pdis.sib.avg <- dat$pdis.sib1
    dat$rev.pdis.sib.avg <- dat$rev.pdis.sib1
  }
  
  l <- lapply(c("avg.maf.IP","n.variants.IP","n.rare.IP","n.causal.IP","pdis.IP","pdis.sib.avg","rev.pdis.sib.avg","num.fam.aff"), function(x) {
    tmp <- data.frame(group=levels(dat$group),
                      measure=x,
                      mean=tapply(dat[,x], dat$group, mean),
                      sd=tapply(dat[,x], dat$group, sd))
    return(tmp)
  })
  
  tmp <- {}
  for(i in 1:length(l)) {
    tmp <- rbind(tmp,l[[i]])
  }
  
  tmp2 <- data.frame(group=levels(dat$group),
                     measure="n",
                     mean=as.numeric(table(dat$group)),
                     sd=NA)
  tmp <- rbind(tmp2,tmp)
  row.names(tmp) <- NULL
  
  
  irep <- unlist(strsplit(file,"individual_data_"))[2]
  tmp$irep <- as.numeric(substr(irep,1,nchar(irep)-4))
  return(dat=tmp)
  
}


##------------------------------------------------------------------------------
## aggregate and print meta output in log txt file
## sim_folder does not need trailing slash
## output_dir does not need trailing slash
##------------------------------------------------------------------------------
#sim_folder <- "/Volumes/JFK_Passport/ordinal_family_history/code/ordinal_flagged_201803121339_post"
#output_dir <- paste0(sim_folder,"/processed_flagged_data")
process_meta_files <- function(sim_folder, output_dir) {
  list_of_folders <- list.files(sim_folder, pattern = "^run_")
  load(paste0(sim_folder,"/parameter_matrix.Rdata"))
  
  ## check
  length(list_of_folders)==nrow(parameter_space)
  
  list_of_folders <- list_of_folders[order(as.numeric(substr(list_of_folders,5,nchar(list_of_folders))))]
  meta <- lapply(list_of_folders, function(x) {
    list_of_files <- list.files(paste0(sim_folder,"/",x), pattern = "^metadata", full.names = TRUE)
    l <- lapply(list_of_files, function(file) {
      mdat <- read.delim(file, stringsAsFactors = F, header = FALSE)
      out <- data.frame(nvariants=as.numeric(mdat$V2[1]),n_rarevariants=as.numeric(mdat$V2[2]),
                        n_causalvariants=as.numeric(mdat$V2[3]),
                        avg_maf_causalvariants=mean(as.numeric(unlist(strsplit(mdat$V2[4],split=",")))),
                        intercept=as.numeric(mdat$V2[7]), beta_geno=as.numeric(mdat$V2[5]),
                        avg_rarebeta=mean(as.numeric(unlist(strsplit(mdat$V2[6],split=",")))))
      return(out)
    })
    mdf <- {}
    for(i in 1:length(l)) {
      mdf <- rbind(mdf, l[[i]])
    }
    mdf <- colMeans(mdf)
    return(mdf)
  })
  
  meta_data <- {}
  for(i in 1:length(meta)) {
    meta_data <- rbind(meta_data, meta[[i]])
  }
  meta_data <- data.frame(meta_data)
  
  sink(paste0(output_dir,"/log.txt"))
  print(colMeans(meta_data))
  sink()
}

##------------------------------------------------------------------------------
## aggregate and save flagged output per group
## sim_folder does not need trailing slash
## output_dir does not need trailing slash
##------------------------------------------------------------------------------
#process_flagged_files(sim_folder="/Volumes/JFK_Passport/ordinal_family_history/code/ordinal_flagged_201805241446_post",
#                      output_dir="/Volumes/JFK_Passport/ordinal_family_history/code/ordinal_flagged_201805241446_post/processed_flagged_data")
process_flagged_files <- function(sim_folder, output_dir) {
  list_of_folders <- list.files(sim_folder, pattern = "^run_")
  load(paste0(sim_folder,"/parameter_matrix.Rdata"))
  
  ## check
  length(list_of_folders)==nrow(parameter_space)
  
  list_of_folders <- list_of_folders[order(as.numeric(substr(list_of_folders,5,nchar(list_of_folders))))]
  
  res <- lapply(list_of_folders, function(x) {
    list_of_files <- list.files(paste0(sim_folder,"/",x), pattern = "^individual_data", full.names = TRUE)
    l <- lapply(list_of_files, read_one_flagged)
    df <- {}
    for(i in 1:length(l)) {
      df <- rbind(df, l[[i]])
    }
    return(df)
  })
  
  names(res) <- list_of_folders
  dat <- {}
  final <- {}
  for(i in 1:length(res)) {
    dat <- res[[i]]
    tmp <- expand.grid(unique(dat$group),unique(dat$measure))
    colnames(tmp) <- c("group","measure")
    for(ii in 1:nrow(tmp)){
      df <- dat[which(dat$group%in%tmp$group[ii] & dat$measure%in%tmp$measure[ii]),]
      tmp[ii,"mean"] <- mean(df$mean)
      tmp[ii,"sd"] <- mean(df$sd)
    }
    tmp <- cbind(tmp,data.frame(parameter_space[i,]))
    final <- rbind(final,tmp)
  }
  
  save(dat, final, res, file=paste0(output_dir,"/flagged.RData"))
}


##------------------------------------------------------------------------------
## make qq plot
##------------------------------------------------------------------------------
myqq <- function(gamut_pvals, unif_pvals, title=NA, truncate = FALSE, ylim=NULL, 
                 thinThreshold=NULL, ci=TRUE, char=20,...) 
{
  pval <- -log10(sort(unif_pvals)) # sort() removes NAs
  n <- length(pval)
  a <- 1:n
  b <- a/n
  x <- -log10(b)
  
  ind <- 1:n
  char <- rep(char,n)
  
  ylm <- ylim
  ylb <- expression(Observed~P~(-log[10]~scale))
  plot(x[ind], pval[ind], type = "n", ylim = ylm, ylab = ylb,
       xlab = expression(Expected~P~(-log[10]~scale)), main=title,cex.main=1,font.lab=3,...)
  # upper and lower have already been subset
  if (ci){
    upper <- qbeta(0.025, a[ind], rev(a)[ind])
    lower <- qbeta(0.975, a[ind], rev(a)[ind])
    polygon(-log10(c(b[ind], rev(b[ind]))), -log10(c(upper, rev(lower))), density=NA, col="gray")
  }
  mx.test <- max(x)
  pval <- -log10(sort(gamut_pvals)) # sort() removes NAs
  n <- length(pval)
  a <- 1:n
  b <- a/n
  x <- -log10(b)
  
  points(x[ind], pval[ind], pch = char, cex=0.7,...)  
  clip(0,mx.test,0,mx.test)
  abline(0,1,col="red",lwd=1.5,xlim=c(0,mx.test))  
  clip(-1,-log10(0.05),-1,-log10(0.05))
  abline(h=-log10(0.05),v=-log10(0.05),col="black",lty=3,lwd=2)
}


##------------------------------------------------------------------------------
## calculate and plot type 1 error
# dir_folder should include trailing slash
##------------------------------------------------------------------------------
#sim_folder <- "/Volumes/JFK_Passport/ordinal_family_history/code/ordinal_error_201806061014_tmp/"
plot_type1_error <- function(dir_folder, run_num=1, ylim=NULL) {
  
  if(is.null(ylim)) {ylim <- c(0,5)}
  
  library(data.table)
  library(reshape2)
  
  #list_of_folders <- list.files(sim_folder, pattern = "^run_")
  load(paste0(dir_folder,"parameter_matrix.Rdata"))
  parameter_space <- parameter_space[run_num,]
  
  #list_of_folders <- list_of_folders[order(as.numeric(substr(list_of_folders,5,nchar(list_of_folders))))]
  sim_folder <- paste0(dir_folder, "run_", run_num)
  list_of_files <- list.files(sim_folder, pattern = "^ordinal_simulation_")
  nsims <- length(list_of_files)
  task_num <- sapply(list_of_files, function(x) substr(x,20,nchar(x)-4))
  
  for(ii in task_num) {
      
      #list_of_files <- list.files(paste0(sim_folder,list_of_folders[i]), pattern = "^ordinal_simulation_")
    
        dat <- read.delim(paste0(sim_folder,"/",paste0("ordinal_simulation_",ii,".txt")),stringsAsFactors=FALSE)
        dat$irep <- paste0("ordinal_simulation_",ii,".txt")
        
        if(exists("aggresults")) {
          aggresults <- rbind(aggresults, dat)
        } else {
          aggresults <- dat
        }
        rm(dat)
  }
  
  save(aggresults,file=paste0(sim_folder,"/","aggresults.Rda"))
  pvals <- aggresults
  pvals <- pvals[,1:3]
  
  #pval_res <- {}  
  #for (i in c(0.01, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4)) {
  #  pval_res <- rbind(pval_res, eval_type1_error(pvals,i))
  #}
  #pval_res <- data.frame(pval_res)
  #pval_res$threshold <- c(0.01, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4)
  
  #pval.idx <- grep("pvalue",colnames(pval_res))
  #melted <- melt(pval_res, id.vars=colnames(pval_res)[-pval.idx], measure.vars=colnames(pval_res)[pval.idx],
  #               variable.name="method", value.name="type_I_error")
  
  
  #library(ggplot2); library(ggthemes)
  #ggplot(data = melted, aes(x = -log10(threshold), y = -log10(type_I_error))) + geom_point() + 
  #  theme_few()  + geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.5) + facet_grid(~method) +
  #  scale_x_continuous("-log10(alpha)") + ylab("-log10(observed type I error)") + 
  #  scale_y_continuous(limits=c(min(-log10(melted$type_I_error)), 6))
  
  library(readr)
  library(plotrix)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(ggpubr)
  
  test <- runif(nsims,0,1)
  
  myqq(gamut_pvals=pvals[,1],unif_pvals=test,title=colnames(pvals)[1],ylim=ylim)
  qq1 <- recordPlot()
  plot.new()
  
  myqq(gamut_pvals=pvals[,2],unif_pvals=test,title=colnames(pvals)[2],ylim=ylim)
  qq2 <- recordPlot()
  plot.new()
  
  myqq(gamut_pvals=pvals[,3],unif_pvals=test,title=colnames(pvals)[3],ylim=ylim)
  qq3 <- recordPlot()
  plot.new()
  
  parameter_space <- cbind(nsims,parameter_space)
  
  return(list(qq1=qq1, qq2=qq2, qq3=qq3, parameter_space=parameter_space))
}
