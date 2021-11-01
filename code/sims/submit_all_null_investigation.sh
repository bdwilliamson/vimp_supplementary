#!/bin/bash

# ------------------------------------------------------------------------------
# submit all null sims
#  1: the simulation name (e.g., "investigate-null-props")
#  2: the learner of interest (e.g., "glm" or "gam")
#  3: the number of total Monte-Carlo reps (e.g., 500)
#  4: the number of replicates per job (e.g., 5)
#  5: whether or not to run cross-fitting (1 is yes, 0 is no)
#  6: which variable to run importance for (1 or 2 [default])
#  7: how many variables to generate (e.g., 2)
#  8: should we generate correlated data? (1 is yes, 0 is no)
#  9: should we use a bootstrap? (1 is yes, 0 is no)
# 10: prefix for i/o files
# 11: should we use the restart queue?
# 12: list of specific jobs to run
# ------------------------------------------------------------------------------
features=(3 2)
estimators=("glm" "gam" "ranger" "SL")
cvs=(0 1)
for feature in ${features[@]}; do
    for cv in ${cvs[@]}; do
        for est in ${estimators[@]}; do
            ./submit_null_investigation.sh "investigate-null-props" $est 1000 5 $cv $feature 4 0 0 "null_${est}_${feature}" 1 ""
        done
    done
done
