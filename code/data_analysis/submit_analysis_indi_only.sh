#!/bin/bash

## submit the vrc01 catnap analysis
ml R/3.5.3-foss-2016b-fh1
ml jbigkit

## takes one command-line arg:
## 1: outcome name (e.g., "cens")

sbatch -M beagle -c10 -p largenode --mem 33G --time=7-0 ./run_analysis_indi_only.sh $1
