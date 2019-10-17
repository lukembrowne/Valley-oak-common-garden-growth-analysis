#!/bin/bash

# Start SNP loop
for n_snps in 25 50 100 200 300 400 500 600 700 800 900 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000 10500 11000 11500 12000 12357;
  do
  
 # If submitting to hoffman
  qsub -cwd -V -N snps_${n_snps} -m bea -M lukembrowne@gmail.com -l h_data=4G,h_rt=24:00:00 -pe shared 8 -o ./logs/snps_${n_snps}.output.txt -e ./logs/snps_${n_snps}.error.txt ./04e-gebv-permutation-test.sh ${n_snps} ./gebv_calculation_notacrossfamily.Rdata

done # End snp loop








