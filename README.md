# Supplementary materials for the general `vimp` paper

This repository contains code to reproduce the analyses in ["A general framework for inference on algorithm-agnostic variable importance"](https://arxiv.org/abs/2004.03683) by Williamson, Gilbert, Simon, and Carone (*Journal of the American Statistical Association (Theory & Methods)*, 2021). All analyses were implemented in the freely available R programming language; specifically, version 4.0.2. All analyses use the R package `vimp` version 2.2.3.

This README file provides an overview of the code available in the repository.  

## Code directory

We have separated our code further into two sub-directories based on the two main objectives of the manuscript:

1. Numerical experiments to evaluate the operating characteristics of our proposed method under varying data-generating mechanisms (`sims`).
2. An analysis of an antibody against HIV-1 infection (`data_analysis`).

All analyses were performed on a Linux cluster using the Slurm batch scheduling system. The head node of the batch scheduler allows the shorthand "ml" in place of "module load". If you use a different batch scheduling system, the individual code files are flagged with the line where you can change batch variables. If you prefer to run the analyses locally, you may -- however, these analyses will then take a large amount of time.

-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/bdwilliamson/vimp_supplementary/issues).
