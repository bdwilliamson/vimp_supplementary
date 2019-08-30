# Running the numerical experiments for the general `vimp` paper

This file describes how to reproduce the simulations in the "A unified approach to nonparametric variable importance assessment" by Williamson, Gilbert, Carone, and Simon. While the code in this file assumes that the user is submitting batch jobs to a high-performance computing (HPC) cluster using the Slurm batch scheduing system, minor edits to these commands allow the use of either local or alternative HPC cluster environments.

Each of the following simulations uses the same R and bash scripts, with varying input arguments. The R functions are:

* `sim_binary_bivariate.R`: provided a simulation name `sim_name` (e.g., "bivariate_loss"), risk type `risk_type` (e.g., "expected_loss" for V-measures), variable importance measure `vimp_measure` (e.g., "deviance"), number of total replicates `nreps_total` (e.g., 1000), number of replicates per HPC job `nreps_per_job` (e.g., 50), number of boostrap replicates `b` (e.g., 1000), and cross-validation flag `cv` (e.g., 1 for using cross-validated VIM estimators), runs the simulation for a specified sample size and variable of interest by replicating `run_sim_binary_bivariate_once.R` a total of `nreps_per_job` times.
* `run_sim_binary_bivariate_once.R`: using arguments specified in `sim_binary_bivariate.R`, runs the simulation a single time; calls `sim_binary_bivariate_data.R` to create the data, and `sim_binary_bivariate_ests.R` to run the estimators.
* `sim_binary_bivariate_data.R`: creates a single dataset of the given sample size from the specified distribution.
* `sim_binary_bivariate_ests.R`: provides the VIM estimators.

The bash scripts are:

*

## Cross-validated plug-in estimators based on flexible techniques

This simulation may be executed from the command line as follows:

```{sh}
./submit_all_sim_binary_bivariate.sh "bivariate_loss" 1000 50 1000 1
./submit_all_sim_binary_bivariate.sh "bivariate_null" 1000 50 1000 1
```

## Plug-in estimators (not cross-validated) based on flexible techniques

This simulation may be executed from the command line as follows:

```{sh}
./submit_all_sim_binary_bivariate.sh "bivariate_loss" 1000 50 1000 0
```

## Effect of not carefully specifying the reward function

This simulation may be executed from the command line as follows:

```{sh}
./submit_all_sim_binary_bivariate.sh "bivariate_naive" 1000 50 1000 0
```