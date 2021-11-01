#!/bin/bash

# submit sims for the simple setting where the learners are either GLM or GAM
# (no Super Learning)

ml fhR/4.0.2-foss-2019b
ml jbigkit

# takes 7 command-line arguments
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
# 11: should we submit to restart? (1 is yes, 0 is no)
# 12: a comma-separated list of jobs that we need to run

# if we haven't passed a list of jobs to run, compute it
if [ "${12}" == "" ]; then
    if [ $7 -ge 5 ]; then
        num_n=2
    else
        if [ $9 -ge 1 ]; then
            num_n=5
        else
            num_n=8
        fi
    fi
    njobs=`expr $3 / $4 \* $num_n`
    arry="1-$njobs"
else
    arry=${12}
fi
io_prefix=${10}
mkdir -p $io_prefix
io_file="$io_prefix/slurm-%A_%a.out"

echo -e '#!/bin/bash\n Rscript investigate_alternative_properties.R' \
        '--sim-name $1 --learner $2 --nreps-total $3 ' \
        '--nreps-per-job $4 --cv $5 --j $6 --p $7' \
        '--corr $8 --bootstrap $9' > investigate_alt.sh
chmod u+x investigate_alt.sh
if [ ${11} -eq 0 ]; then
    sbatch -A gilbert_p --time=1-0 --array=$arry -e $io_file -o $io_file ./investigate_alt.sh $1 $2 $3 $4 $5 $6 $7 $8 $9
else
    sbatch --qos=restart-new --partition=restart-new --time=1-0 --array=$arry -e $io_file -o $io_file ./investigate_alt.sh $1 $2 $3 $4 $5 $6 $7 $8 $9
fi
rm investigate_alt.sh
