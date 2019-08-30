#!/bin/bash

## takes 6 command-line args:
## 1: the simulation name (e.g., "bivariate")
## 2: the risk type (e.g., "expected_loss")
## 3: the measure of importance (e.g., "deviance")
## 4: the number of total Monte-Carlo reps (e.g., 1000)
## 5: the number of replicates per batch job (e.g., 50)
## 6: the number of bootstrap replicates (e.g., 1000)

Rscript sim_binary_bivariate.R --sim-name $1 --risk-type $2 --vimp-measure $3 --nreps-total $4 --nreps-per-job $5 --b $6