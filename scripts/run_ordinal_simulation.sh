#!/bin/sh 
# Tell the SGE that this is an array job, with 'tasks' to be numbered 1 to 10000 
#$ -t 1-10000 
# When a single command in the array job is sent to a compute node, 
# its task number is stored in the variable SGE_TASK_ID, 
# so we can use the value of that variable to get the results we want. 
 
module load R 
module load cosi/1.2 
 
PRJDIR="${HOME}/mpe_projects/ordinal_gamut/ordinal_error_201806181156/run_1" 
 
if [ -e /scratch/$JOB_ID ]; then 
     echo "Error. Job folder (/scratch/$JOB_ID) exists on node `hostname`" 
     exit 
else 
     mkdir /scratch/$JOB_ID 
     SCRATCHDIR="/scratch/$JOB_ID" 
fi 
 
cp ${PRJDIR}/* ${SCRATCHDIR} 
cd ${SCRATCHDIR} 
 
Rscript ${SCRATCHDIR}/03_simulation_ordinal_one.r $SGE_TASK_ID 
 
/bin/rm -f ${SCRATCHDIR}/autosomes_deCODE.distr 
/bin/rm -f ${SCRATCHDIR}/model.test 
/bin/rm -f ${SCRATCHDIR}/out.* 
/bin/rm -f ${SCRATCHDIR}/params 
/bin/rm -f ${SCRATCHDIR}/recParams 
/bin/rm -f ${SCRATCHDIR}/logfile.txt 
/bin/rm -f ${SCRATCHDIR}/run.pl 
/bin/rm -f ${SCRATCHDIR}/*.r 
/bin/rm -f ${SCRATCHDIR}/*.sh 
/bin/rm -f ${SCRATCHDIR}/*.R* 
 
rsync -av ${SCRATCHDIR}/ ${PRJDIR} 
 
/bin/rm -fr ${SCRATCHDIR} 
 
module unload R 
module unload cosi/1.2
