#!/bin/bash

# Command to submit job
# qsub -cwd -l h_data=8G,h_rt=05:00:00 -t 2-20:4 run_simulation_total_analysis.sh

# qsub -cwd -l h_data=8G,h_rt=05:00:00 run_simulation_total_analysis.sh

## Script to generate simulated datasets, 
## analyze them with GAM and BLUP framework, and compare, all on the cluster


# Load modules
source /u/local/Modules/default/init/modules.sh ## Need .sh ending since this is a shell script
module load R/3.5.1
module li
 
. /u/local/etc/profile.d/sge.sh # To avoid qsub command not found error

n_sites=$1
qtl_ve=$2
rep=$3

echo $n_sites
echo $qtl_ve
echo $rep

# Set number of sites to simulate
#n_sites=${n_sites}

## If submitting as a job array
#n_sites=$SGE_TASK_ID

# Run suffix
suffix="${n_sites}_sites_${qtl_ve}_rep${rep}"


# Make and change directory
mkdir -p sim_output_${suffix}
cp ./08_simulating_data.R ./sim_output_${suffix}
cp ./08_analyzing_simulated_data.R ./sim_output_${suffix}
cp ./run_gams_sim.R ./sim_output_${suffix}
cp ./run_gams_sim.R.cmd ./sim_output_${suffix}
cd ./sim_output_${suffix} # Change directory

mkdir -p figures
mkdir -p bglr

## First, simulate dataset
 R CMD BATCH --vanilla "--args $n_sites $qtl_ve $rep" 08_simulating_data.R

 # Rename Rplots.pdf file
 mv Rplots.pdf Rplots_simulation.pdf

## Need to modify command script to add in specific variables before submitting as job array
sed -i "s/#set run_label=\"XXX\"/set run_label=\"${n_sites}sites\"/g" run_gams_sim.R.cmd
sed -i "s/#set data_file=\"XXX\"/set data_file=\"gam_cluster_sim_${n_sites}sites_workspace.Rdata\"/g" run_gams_sim.R.cmd


# Submit job array
  qsub -sync y -N r${rep}_${n_sites}_${qtl_ve} run_gams_sim.R.cmd
  wait

 ## Get path to GAM output - should only get gam results folder
	 gam_output="$(find . -name "run*" -type d)"   
	 data_file=gam_cluster_sim_${n_sites}sites_workspace.Rdata 

## Run code to analyze simulated data
R CMD BATCH --vanilla "--args $gam_output $data_file $n_sites" 08_analyzing_simulated_data.R

echo "All done."









