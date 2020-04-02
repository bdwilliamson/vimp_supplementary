#!/bin/bash

Rscript sim_no_sample_splitting.R --sim-name $1 --risk-type "expected_loss" --vimp-measure "accuracy" "auc" --nreps-total $2 --nreps-per-job $3 --b $4 --cv $5 --delta $6
