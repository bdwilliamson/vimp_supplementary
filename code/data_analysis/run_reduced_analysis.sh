#!/bin/bash

## run reduced analysis of only subtype
Rscript run_super_learners.R --outcome $1 -redgrps -rungeog -runindi
