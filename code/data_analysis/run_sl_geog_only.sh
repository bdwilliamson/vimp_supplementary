#!/bin/bash

## run geog only analysis
Rscript run_super_learners.R --outcome $1 -runfull -rungrps -runindi
