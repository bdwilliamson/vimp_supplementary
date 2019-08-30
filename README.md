# Supplementary materials for the general `vimp` paper

This repository contains code to reproduce the analyses in "A unified approach to nonparametric variable importance assessment" by Williamson, Gilbert, Carone, and Simon. All analyses were implemented in the freely available R programming language; specifically, version 3.5.3. All analyses use the R package `vimp`, version 2.0.0.

This README file provides an overview of the code available in the repository.  

## Code directory

We have separated our code further into two sub-directories based on the two main objectives of the manuscript:

1. Numerical experiments to evaluate the operating characteristics of our proposed method under varying data-generating mechanisms (`sims`).
2. An analysis of an antibody against HIV-1 infection (`data_analysis`).

Within each sub-directory, we further subdivide the code into an R directory (hosting all of the R code for the analysis) and a Shell directory (hosting all of the code for batch submission to a high-performance cluster computing environment). All analyses were performed on a Linux cluster using the Slurm batch scheduling system. If you use a difference batch scheduling system, the individual code files are flagged with the line where you can change batch variables. If you prefer to run the analyses locally, you may -- however, these analyses will then take a large amount of time.

-----

## Issues

If you encounter any bugs or have any specific questions about the analysis, please
[file an issue](https://github.com/bdwilliamson/vimp_supplementary/issues).
