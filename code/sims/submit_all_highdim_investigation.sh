#!/bin/bash

# ------------------------------------------------------------------------------
# submit all null-hypothesis/alt-hypothesis sims in higher-dimensions
#  1: the simulation name (e.g., "investigate-alt-props")
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
# dimensions=(50 100 200)
# features=("1" "2" "3" "4" "13" "24")
# corrs=(0 1)
# for dim in ${dimensions[@]}; do
#     for corr in ${corrs[@]}; do
#         for feature in "${features[@]}"; do
#             ./submit_alt_investigation.sh "investigate-highdim-props" "SL" 1000 2 1 $feature $dim $corr 0 "highdim_${dim}_${feature}" 1 ""
#         done
#     done
# done
dimensions=(50 100 200)
features=("1" "2" "3" "4" "13" "24")
corrs=(0 1)
for dim in ${dimensions[@]}; do
    for corr in ${corrs[@]}; do
        for feature in "${features[@]}"; do
            ./submit_null_investigation.sh "investigate-highdim-props" "SL" 1000 2 1 $feature $dim $corr 0 "highdim_${dim}_${feature}" 1 ""
        done
    done
done
