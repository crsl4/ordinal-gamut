load("parameter_matrix.Rdata")

files_to_copy <- list.files(".")

scenarios <- paste0("run_",1:nrow(parameter_space))

sapply(scenarios,dir.create)

write_run_dirs <- function(x) {
  file.copy(from=files_to_copy, to=paste0("./",scenarios[x]),overwrite=TRUE)
  params <- parameter_space[x,]
  assign("seqlength",as.double(params[1]))
  assign("nhap_eur",as.double(params[2]))
  assign("nhap_afr",as.double(params[3]))
  assign("n_per_group",as.double(params[4]))
  assign("Nsibsize",as.double(params[5]))
  assign("minaff",as.double(params[6]))
  assign("include_parents",as.double(params[7]))
  assign("randomize_nsibs",as.double(params[8]))
  assign("rare_maf",as.double(params[9]))
  assign("p_causal",as.double(params[10]))
  assign("p_adverse",as.double(params[11]))
  assign("effectsize",as.double(params[12]))
  assign("prev_target",as.double(params[13]))
  assign("lambda",as.double(params[14]))
  assign("p_yri",as.double(params[15]))
  assign("beta_yri",as.double(params[16]))
  assign("cconfounder_family",as.double(params[17]))
  assign("beta_cconfounder",as.double(params[18]))
  assign("residual_format",as.character(params[19]))
  
  rm(params,envir=environment())
  save(list=ls(all.names=TRUE),file=paste0("./",scenarios[x],"/arguments.Rda"),envir=environment())
  
  text = paste("\n seqlength =",seqlength,
               "\n nhap_eur =",nhap_eur,
               "\n nhap_afr =",nhap_afr,
               "\n n_per_group =",n_per_group,
               "\n Nsibsize =",Nsibsize,
               "\n minaff =",minaff,
               "\n include_parents =",include_parents,
               "\n randomize_nsibs =",randomize_nsibs,
               "\n rare_maf =",rare_maf,
               "\n p_causal =",p_causal,
               "\n p_adverse =",p_adverse,
               "\n effectsize =",effectsize,
               "\n prev_target =",prev_target,
               "\n lambda =",lambda,
               "\n p_yri =",p_yri,
               "\n beta_yri =",beta_yri,
               "\n cconfounder_family =",cconfounder_family,
               "\n beta_cconfounder =",beta_cconfounder,
               "\n residual_format =",residual_format
  )
  write(text,paste0("./",scenarios[x],"/logfile.txt"))
}

sapply(1:length(scenarios),write_run_dirs)

write_submit_scripts <- function(x, nsims, path_to_parent_dir) {
  text <- paste0(
    "#!/bin/sh \n",
    "# Tell the SGE that this is an array job, with 'tasks' to be numbered 1 to 10000 \n",
    paste0("#$ -t 1-",nsims," \n"),
    "# When a single command in the array job is sent to a compute node, \n",
    "# its task number is stored in the variable SGE_TASK_ID, \n",
    "# so we can use the value of that variable to get the results we want. \n \n",
    
    "module load R \n",
    "module load cosi/1.2 \n \n",
    
    paste0('PRJDIR="${HOME}/',path_to_parent_dir,"/",scenarios[x],'" \n \n'),
    
    
    "if [ -e /bin/mktemp ]; then \n",
    "TMPDIR=`/bin/mktemp -d /scratch/XXXXXX` \n",
    "elif [ -e /usr/bin/mktemp ]; then \n",
    "TMPDIR=`/usr/bin/mktemp -d /scratch/XXXXXX` \n",
    "else \n",
      "echo 'Error. Cannot find program to create tmp directory' \n",
    "exit \n",
    "fi \n \n",
    
    "cp ${PRJDIR}/* ${TMPDIR} \n",
    "cd ${TMPDIR} \n \n",
    
    "Rscript ${TMPDIR}/03_simulation_ordinal_one.r $SGE_TASK_ID \n \n",
    
    "/bin/rm -f ${TMPDIR}/autosomes_deCODE.distr \n",
    "/bin/rm -f ${TMPDIR}/model.test \n",
    "/bin/rm -f ${TMPDIR}/out.* \n",
    "/bin/rm -f ${TMPDIR}/params \n",
    "/bin/rm -f ${TMPDIR}/recParams \n",
    "/bin/rm -f ${TMPDIR}/logfile.txt \n",
    "/bin/rm -f ${TMPDIR}/run.pl \n",
    "/bin/rm -f ${TMPDIR}/*.r \n",
    "/bin/rm -f ${TMPDIR}/*.sh \n",
    "/bin/rm -f ${TMPDIR}/*.Rda \n \n",
    
    "rsync -av ${TMPDIR}/ ${PRJDIR} \n \n",
    
    "/bin/rm -fr ${TMPDIR} \n \n",
    
    "module unload R \n",
    "module unload cosi/1.2"
  )
  
  write(text, paste0("./",scenarios[x],"/run_ordinal_simulation.sh"))
}

##### entire desired path here ####
#sapply(1:length(scenarios), write_submit_scripts, nsims=10000,
#       path_to_parent_dir="mpe_projects/ordinal_gamut/ordinal_error_201806061014")


