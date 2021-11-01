#!/bin/bash

## run vrc01 catnap analysis

## ----------------------------------------------------
## create the dataset
## ----------------------------------------------------
Rscript get_catnap_data.R

## ----------------------------------------------------
## run the super learners
## ----------------------------------------------------
Rscript run_super_learners.R --outcome $1 --ind-type $2

## ----------------------------------------------------
## get variable importance
## ----------------------------------------------------
Rscript get_vimp.R --outcome $1 --vimp-measure "r_squared"

Rscript get_vimp.R --outcome $1 --vimp-measure "auc"

Rscript get_vimp.R --outcome $1 --vimp-measure "accuracy"
