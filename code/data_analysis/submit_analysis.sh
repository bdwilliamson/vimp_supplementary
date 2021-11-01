#!/bin/bash

## submit the vrc01 catnap analysis
ml fhR/4.0.2-foss-2019b
ml jbigkit

## takes command-line args:
## 1: outcome name (e.g., "cens")
## 2: individual importance type ("sitewise" [AA sites] or "residuewise" [residues at AA sites])

sbatch -A gilbert_p -c10 --mem 33G --time=7-0 ./run_analysis.sh $1 $2
