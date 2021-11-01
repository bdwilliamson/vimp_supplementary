# Running the data analysis for the general `vimp` paper

This file describes how to reproduce the data analysis in ["A general framework for inference on algorithm-agnostic variable importance"](https://arxiv.org/abs/2004.03683) by Williamson, Gilbert, Simon, and Carone (*Journal of the American Statistical Association (Theory & Methods)*, 2021). All analyses were implemented in the freely available R programming language; specifically, version 4.0.2 All analyses use the R package `vimp` version 2.2.3.

The R scripts contain the code necessary to compute estimates:
* `get_catnap_data.R`: get the CATNAP data used in Magaret, Benkeser, Williamson et al. (*PLoS Computational Biology*, 2019)
* `run_super_learners.R`: run the specified Super Learner estimator on the CATNAP data, for the specified outcome
* `get_vimp.R`: obtain estimates of variable importance, defined using R-squared, AUC, and classification accuracy
* `render_plots.R`: generate plots based on a variable importance analysis
* `get_variable_groups.R`: define groups of variables (amino acid sequence and viral geometry features)
* `super_learner_libraries.R`: functions to define the Super Learner library
* `utils.R`: various utility functions

The bash scripts run the analysis:
* `run_analysis.sh`: run the full analysis for the specified outcome
* `submit_analysis.sh`: submit batch jobs for full analysis
* `render_plots.sh`: render plots once we have the required output

The analyses presented in the manuscript can be replicated using the following code (note that this may take some time to run, especially if you aren't using a Slurm scheduler):
```{sh}
chmod u+x *.sh
# the main analysis
./run_analysis.sh "sens50" "sitewise"
# harmonized with Magaret et al. (2019)
./run_analysis.sh "cens" "sitewise"
# render all of the plots
./render_plots.sh
```
