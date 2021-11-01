# Running the numerical experiments for the general `vimp` paper

This file describes how to reproduce the simulations in ["A general framework for inference on algorithm-agnostic variable importance"](https://arxiv.org/abs/2004.03683) by Williamson, Gilbert, Simon, and Carone (*Journal of the American Statistical Association (Theory & Methods)*, 2021). All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. All analyses use the R package `vimp` version 2.2.3.

The numerical experiments consist of four sections: first, we consider the properties of our proposal under the alternative hypothesis (Scenario 1, with two important features); next, we consider the properties of our proposal under the null hypothesis (Scenario 2, with four features, two of which are important and two of which are unimportant). We next consider a bootstrap for interval estimation under the alternative hypothesis. Finally, we investigate the properties of our proposal with a larger number of features (both with and without correlation).

To obtain the true values of accuracy and AUC under our simulation model, we used the file `compute_truths_probit.R`.

## Properties of our proposal under the alternative hypothesis

The following code will replicate the results in Section 5.2 of the main manuscript (Figure 2) and Section 4.2 of the Supplementary Material (Figures S1--S4).

The simulation uses the following files:
* `submit_all_alt_investigation.sh`: Submit all simulations for this section.
* `submit_alt_investigation.sh`: Batch submission of a group of jobs to a Slurm scheduler.
* `investigate_alternative_properties.R`: the main R script for this simulation. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `investigate_alternative_once.R`: Runs the simulation a single time for a specified set of parameters.
* `gen_data.R`: Generate a dataset.
* `run_ests.R`: Run a single estimator (SuperLearner, random forests, logistic regression, GAMs, etc.) with or without cross-fitting and sample-splitting
* `utils.R`: Utility functions.

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_all_alt_investigation.sh
```
If you aren't running on a Slurm scheduler, make sure to edit the appropriate lines (flagged in each file). You can run this code locally, but it will take some time.

## Properties of our proposal under the null hypothesis

The following code will replicate the results in Section 5.2 of the main manuscript (Figure 3) and Section 4.3 of the Supplementary Material (Figures S5--S8).

The simulation uses the following files:
* `submit_all_null_investigation.sh`: Submit all simulations for this section.
* `submit_null_investigation.sh`: Batch submission of a group of jobs to a Slurm scheduler.
* `investigate_null_properties.R`: the main R script for this simulation. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `investigate_null_once.R`: Runs the simulation a single time for a specified set of parameters.
* `gen_data.R`: Generate a dataset.
* `run_ests.R`: Run a single estimator (SuperLearner, random forests, logistic regression, GAMs, etc.) with or without cross-fitting and sample-splitting
* `utils.R`: Utility functions.

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_all_null_investigation.sh
```
If you aren't running on a Slurm scheduler, make sure to edit the appropriate lines (flagged in each file). You can run this code locally, but it will take some time.

## Using the bootstrap for interval estimation

The following code will replicate the results in Section 4.4 of the Supplementary Material (Figures S9--S12).

The simulation uses the following files:
* `submit_all_boot_investigation.sh`: Submit all simulations for this section.
* `submit_alt_investigation.sh`: Batch submission of a group of jobs to a Slurm scheduler.
* `investigate_alternative_properties.R`: the main R script for this simulation. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `investigate_alternative_once.R`: Runs the simulation a single time for a specified set of parameters.
* `gen_data.R`: Generate a dataset.
* `run_ests.R`: Run a single estimator (SuperLearner, random forests, logistic regression, GAMs, etc.) with or without cross-fitting and sample-splitting
* `utils.R`: Utility functions.

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_all_boot_investigation.sh
```
If you aren't running on a Slurm scheduler, make sure to edit the appropriate lines (flagged in each file). You can run this code locally, but it will take some time.

## Higher dimensions and correlated features

The following code will replicate the results in Section 4.5 of the Supplementary Material (Figures S13--S16).

The simulation uses the following files:
* `submit_all_highdim_investigation.sh`: Submit all simulations for this section.
* `submit_null_investigation.sh`: Batch submission of a group of jobs to a Slurm scheduler.
* `investigate_null_properties.R`: the main R script for this simulation. Runs the simulation `nreps_per_job` times for a specified set of parameters.
* `investigate_null_once.R`: Runs the simulation a single time for a specified set of parameters.
* `gen_data.R`: Generate a dataset.
* `run_ests.R`: Run a single estimator (SuperLearner, random forests, logistic regression, GAMs, etc.) with or without cross-fitting and sample-splitting
* `utils.R`: Utility functions.

Running the following code will submit all of the simulations to a Slurm scheduler:
```{bash}
chmod u+x *.sh
./submit_all_highdim_investigation.sh
```
If you aren't running on a Slurm scheduler, make sure to edit the appropriate lines (flagged in each file). You can run this code locally, but it will take some time.

## Compiling results and making figures

Once you have results from all of the simulations, you can run the following code to create all of the figures:
```{bash}
chmod u+x *.sh
./create_all_plots.sh
```
This code calls, in turn:
* `load_sim_alt_lowdim.R`: load simulations from the alternative hypothesis simulation and create plots
* `load_sim_bootstrap.R`: load simulations from the bootstrap simulation and create plots
* `load_sim_null.R`: load simulations from the null hypothesis simulation and create plots
* `load_sim_highdim.R`: load simulations from the higher-dimensional simulation and create plots
* `create_combined_plots.R`: combine plots together for the supplement
