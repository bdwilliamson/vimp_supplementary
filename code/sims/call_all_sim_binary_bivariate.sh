#!/bin/bash

## takes 6 command-line args:
## 1: the simulation name (e.g., "bivariate")
## 2: the number of total Monte-Carlo reps (e.g., 1000)
## 3: the number of replicates per batch job (e.g., 50)
## 4: the number of bootstrap replicates (e.g., 1000)

Rscript sim_binary_bivariate.R --sim-name $1 --risk-type "expected_loss" --vimp-measure "deviance" "accuracy" "auc" --nreps-total $2 --nreps-per-job $3 --b $4 --cv $5
