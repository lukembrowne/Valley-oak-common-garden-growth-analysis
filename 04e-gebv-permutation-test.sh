#!/bin/bash

 source /u/local/Modules/default/init/modules.sh ## Need .sh ending since this is a shell script
 module load R/3.5.1

Rscript ./04e-gebv-permutation-test.R $1 $2