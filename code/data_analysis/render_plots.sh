#!/bin/bash

## generate plots from VRC01 vimp analysis; note .05 / 13 groups = this threshold
Rscript render_plots.R --outcome "cens" "sens50" --ind-type "sitewise" --threshold 0.003846154
