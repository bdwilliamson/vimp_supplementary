# Running the data analysis for the general `vimp` paper

This file describes how to reproduce the data analysis in the "A unified approach to nonparametric variable importance assessment" by Williamson, Gilbert, Carone, and Simon. The code in this file uses the prediction functions estimated as a result of "Prediction of VRC01 neutralization sensitivity by HIV-1 gp160 sequence features" by Magaret, Benkeser, Williamson et al. (2019). Please run all code from [this GitHub repository](https://github.com/benkeser/vrc01) prior to running any analyses for the current paper.

Running `vrc01_vimp_analysis_compute_group_vimp.R` computes variable importance estimates for the groups defined by Magaret et al. (2019) based on the fitted values resulting from the analysis of Magaret et al. (2019), but using the deviance, accuracy, and AUC.

Running `vrc01_vimp_analysis_compute_group_vimp.R` computes variable importance estimates for the individual features defined by Magaret et al. (2019) based on the fitted values resulting from the analysis of Magaret et al. (2019), but using the deviance, accuracy, and AUC.


Running `vrc01_vimp_analysis.R` loads in the results of the previous two function calls and creates plots and tables of results.

These three functions rely on `vrc01_vimp_analysis_helpers.R` for utility functions.

The final plots for the paper are created using `vrc01_vimp_analysis_final_plots.R`.
