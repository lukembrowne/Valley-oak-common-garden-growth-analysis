#!/bin/bash
for rep in 13
do
for n_sites in 10 20
	do
	for qtl_ve in 0.01 0.99
		do
		# Submit job to cluster	
		qsub -cwd -N r${rep}_${n_sites}_${qtl_ve} -l h_data=16G,h_rt=24:00:00 run_simulation_total_analysis.sh $n_sites $qtl_ve $rep

	done # End qtl ve loop
	done # End site loop
done # End Rep loop

