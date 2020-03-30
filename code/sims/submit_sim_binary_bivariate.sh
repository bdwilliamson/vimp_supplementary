#!/bin/bash

ml R/3.5.3-foss-2016b-fh1
ml jbigkit

## takes 6 command-line args:
## 1: the simulation name (e.g., "bivariate")
## 2: the risk type (e.g., "expected_loss")
## 3: the measure of importance (e.g., "deviance")
## 4: the number of total Monte-Carlo reps (e.g., 1000)
## 5: the number of replicates per batch job (e.g., 50)
## 6: the number of bootstrap replicates (e.g., 1000)

## num_n_j is the number of unique n, j combinations (10 ns * 2 js for small n simulation)
num_n_j=20
njobs=`expr $4 / $5 \* $num_n_j`
echo $njobs

sbatch --array=1-$njobs ./call_sim_binary_bivariate.sh $1 $2 $3 $4 $5 $6
