#!/bin/csh -f
#  run_gams.R.cmd
#
#  UGE job for run_gams.R built Tue May 22 10:49:43 PDT 2018
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/flashscratch/flashscratch2/l/lukembro/qlobata_growth/run_gams.R.joblog.$JOB_ID.$TASK_ID
#$ -o /u/flashscratch/flashscratch2/l/lukembro/qlobata_growth/logs/run_gams.R.joblog.$JOB_ID.$TASK_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/flashscratch/flashscratch2/l/lukembro/qlobata_growth/run_gams.R
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Threaded:     8-way threaded
#  Resources requested
#$ -l h_data=4096M,h_rt=01:00:00
# #
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M lukembro@mail
#  Notify at beginning and end of job
#$ -m a
#  Job is not rerunable
#$ -r n
#
#  Job array indexes
#$ -t 1-16979:100
#
# 16979 Snps after iterative filtering

# Set variable that will be passed to R script
# Interval needs to match what is set in Job array indexes

set run_label="rgr_fREML_discrete_factor"
set interval=100
set climate_var_dif='tmax_sum_dif'



# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "jobarray threads"
  set qqmtasks  = 8
  set qqidir    = /u/flashscratch/flashscratch2/l/lukembro/qlobata_growth
  set qqjob     = run_gams.R
  set qqodir    = /u/flashscratch/flashscratch2/l/lukembro/qlobata_growth
  cd     /u/flashscratch/flashscratch2/l/lukembro/qlobata_growth
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for run_gams.R built Tue May 22 10:49:43 PDT 2018"
  echo ""
  echo "  run_gams.R directory:"
  echo "    "/u/flashscratch/flashscratch2/l/lukembro/qlobata_growth
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
  echo "  run_gams.R 8-way threaded job configuration:"
#
  echo ""
  echo "Task $SGE_TASK_ID for run_gams.R started on:   "` hostname -s `
  echo "Task $SGE_TASK_ID for run_gams.R started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load R/3.4.2
  module li
#

# Make directories and copy over R file
mkdir run_${JOB_ID}_${climate_var_dif}_{$run_label}
cp ./run_gams.R ./run_${JOB_ID}_${climate_var_dif}_{$run_label}
cp ./run_gams.R.cmd ./run_${JOB_ID}_${climate_var_dif}_{$run_label}
cp ./gam_cluster_*.Rdata ./run_${JOB_ID}_${climate_var_dif}_{$run_label}
cd ./run_${JOB_ID}_${climate_var_dif}_{$run_label}
mkdir R_output
mkdir output
mkdir gam_mods_data
mkdir model_summaries
mkdir model_plots

#
  # if more than 10 jobtasks, send only "a" mail.  until...
  # ... the last but 1 task so it applies to last task. then send or not per user preference.
  if( 2 > 10 && $SGE_TASK_ID == 5 ) then
    /u/local/bin/qalter -m a $JOB_ID
    sleep 2
  endif
#
  echo run_gams.R "" \>\& run_gams.R.output.$JOB_ID.$SGE_TASK_ID
  echo ""
  setenv OMP_NUM_THREADS 8
#
  /usr/bin/time R CMD BATCH --vanilla "--args $SGE_TASK_ID $interval $climate_var_dif" run_gams.R ./R_output/run_gams.$JOB_ID.$SGE_TASK_ID.out >& ./output/run_gams.R.output.$JOB_ID.$SGE_TASK_ID
#
  echo ""
  echo "Task $SGE_TASK_ID for run_gams.R finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/flashscratch/flashscratch2/l/lukembro/qlobata_growth/logs/run_gams.R.joblog.$JOB_ID.$SGE_TASK_ID --------" >> /u/local/apps/queue.logs/jobarray.log.multithread
  if (`wc -l ./logs/run_gams.R.joblog.$JOB_ID.$SGE_TASK_ID  | awk '{print $1}'` >= 1000) then
        head -50 ./logs/run_gams.R.joblog.$JOB_ID.$SGE_TASK_ID >> /u/local/apps/queue.logs/jobarray.log.multithread
        echo " "  >> /u/local/apps/queue.logs/jobarray.log.multithread
        tail -10 ./logs/run_gams.R.joblog.$JOB_ID.$SGE_TASK_ID >> /u/local/apps/queue.logs/jobarray.log.multithread
  else
        cat ./logs/run_gams.R.joblog.$JOB_ID.$SGE_TASK_ID >> /u/local/apps/queue.logs/jobarray.log.multithread
  endif
  exit (0)
