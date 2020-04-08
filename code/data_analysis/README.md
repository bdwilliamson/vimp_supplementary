# Running the data analysis for the general `vimp` paper

This file describes how to reproduce the data analysis in "A unified approach for inference on algorithm-agnostic variable importance" by Williamson, Gilbert, Simon, and Carone (*arXiv*, 2020). The code in this file uses the prediction functions estimated as a result of "Prediction of VRC01 neutralization sensitivity by HIV-1 gp160 sequence features" by Magaret, Benkeser, Williamson et al. (2019). All analyses are run using R version 3.5.3 and the `vimp` package version 2.0.1.

The R scripts contain the code necessary to compute estimates:
* `get_catnap_data.R`: get the CATNAP data used in Magaret, Benkeser, Williamson et al. (*PLoS Computational Biology*, 2019)
* `run_super_learners.R`: run the specified Super Learner estimator on the CATNAP data, for the specified outcome
* `get_vimp.R`: obtain estimates of variable importance, defined using R-squared, AUC, and classification accuracy
* `render_plots.R`: generate plots based on a variable importance analysis
* `super_learner_libraries.R`: functions to define the Super Learner library
* `utils.R`: various utility functions

The bash scripts run the analysis:
* `run_analysis.sh`: run the full analysis for the specified outcome
* `run_analysis_indi_only.sh`: run the analysis only for individual variable importance
* `run_reduced_analysis.sh`: run for a specified subgroup of groups
* `run_sl_geog_only.sh`: run for only geographic confounders
* `submit_analysis_indi_only.sh`: submit batch jobs for individual VIM analysis
* `submit_analysis.sh`: submit batch jobs for full analysis
* `render_plots.sh`: render plots once we have the required output

The analyses presented in the manuscript can be replicated using the following code (note that this may take some time to run):
```{sh}
./run_analysis.sh "cens"
./render_plots.sh
```
