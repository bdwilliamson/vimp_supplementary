#!/bin/bash

ml R/3.5.3-foss-2016b-fh1
ml jbigkit

## takes 6 command-line args:
## 1: the simulation name (e.g., "bivariate")
## 2: the number of total Monte-Carlo reps (e.g., 1000)
## 3: the number of replicates per batch job (e.g., 50)
## 4: the number of bootstrap replicates (e.g., 1000)
## 5: whether or not to run cv (1 is run, 0 is not)

## num_n_j is the number of unique n, j combinations (16 ns * 2 js for small n simulation)
num_n_j=20
njobs=`expr $2 / $3 \* $num_n_j`
echo $njobs

sbatch --array=1-$njobs ./call_all_sim_binary_bivariate.sh $1 $2 $3 $4 $5
