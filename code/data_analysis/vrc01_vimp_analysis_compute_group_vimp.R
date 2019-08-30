#!/usr/local/bin/Rscript
##########################################################################
##
## FILE: vrc01_vimp_analysis_compute_group_vimp.R
##
## PURPOSE: variable importance analysis of binary ic50 endpoints
##          from VRC01 CATNAP analysis, using deviance, accuracy, auc
##
##          this file computes variable importance estimates based on
##          estimated regression functions, and saves the main output
##          FOR GROUPS ONLY
##
##########################################################################

##------------------------------------------------------------------------
## set up
##------------------------------------------------------------------------
## load required functions and packages
# only run this if necessary to update package
# install.packages("~/Projects/UW/vimp_1.3.0.tar.gz", repos = NULL, type = "source")
library("vimp")
library("argparse")
code_dir <- "../code/"
plots_dir <- "../code/plots/"
output_dir <- "../code/vrc01_analysis_output/"
source(paste0(code_dir, "vrc01_vimp_analysis_helpers.R"))

parser <- ArgumentParser()
parser$add_argument("--vimp-measure", default = "deviance", help = "the variable importance measure of interest")
parser$add_argument("--level", type = "double", default = 0.95, help = "nominal level for confidence intervals")
args <- parser$parse_args()

## set directories, dataset version, etc.
vrc01_results_dir <- "~/Projects/VIDD/hvtn/vrc01_genotypic_resistance_score/results/"
vrc01_code_dir <- "~/Projects/VIDD/hvtn/vrc01_genotypic_resistance_score/code/"
vrsn <- "12"
nms_pub <- list(c(expression("Sensitive/Resistant Only; deviance"),
             expression(paste(IC[50], " Censored; deviance", sep = ""))),
             c(expression("Sensitive/Resistant Only; accuracy"),
             expression(paste(IC[50], " Censored; accuracy", sep = ""))),
             c(expression("Sensitive/Resistant Only; AUC"),
             expression(paste(IC[50], " Censored; AUC", sep = ""))))
run <- "standard"

nms_group <- c("VRC01 contact sites (Set 1)", "CD4 binding sites (Set 2)", "ESA (Set 3)",
                "Glycosylation sites (Set 4)", "Covarying sites (Set 5)", "PNG sites (Set 6)",
                "gp41 sites (Set 7)",
                "N-linked glycosylation (Set 8)", "Subtype (Set 9)", "Sequons (Set 10)",
                "Viral geometry (Set 11)", "Cysteines (Set 12)", "Steric bulk (Set 13)")
num_grp <- 13

##------------------------------------------------------------------------
## get estimates, CIs, hypothesis test reject/not
##------------------------------------------------------------------------
## -------------------------------------
## fitted values
## -------------------------------------
## load in predicted values from regression on all features
full_mat <- expand.grid(o = c("sens.resis", "cens"),
                          s = 1:2, stringsAsFactors = FALSE)
full_list <- as.list(paste0(vrc01_results_dir, run, "/preds_", full_mat$o, "_set",
                              full_mat$s, "_v", vrsn, "_", run,".rds"))

## read in the data
tmp <- lapply(full_list, readRDS)
fulls <- lapply(lapply(tmp, function(x) x$pred), rev) # need to reverse the inner order due to cvma fitting; see get_preds.R comments for more detail
folds <- lapply(tmp, function(x) x$folds)

## set up the groups
grp_mat <- expand.grid(o = c("sens.resis", "cens"),
                     f = 1:num_grp, s = 1:2, stringsAsFactors = FALSE)
grp_mat_chr <- grp_mat
grp_mat_chr$f <- as.character(grp_mat$f)
grp_mat_chr$s <- as.character(grp_mat$s)

## names of fitted values
red_grp_list <- as.list(paste0(vrc01_results_dir, run, "/reduced_preds_",
                             grp_mat$o, "_dataset_", grp_mat$s, "_group_",
                             grp_mat$f, "_", run,".rds"))

red_grps_rev <- lapply(red_grp_list, readRDS)
grp_nms <- apply(grp_mat_chr, 1, function(x) paste0("grp_", paste(x, collapse = "_")))

## reverse them all; the 10th-fold-fit is currently in the first position for each list
red_grps <- lapply(red_grps_rev, rev)

## create lists with the appropriate full fits, same length as both group and individual ones
fulls_grp <- c(rep(fulls[1:length(unique(grp_mat$o))], dim(grp_mat)[1]/length(unique(grp_mat$o))/2),
                 rep(fulls[(length(unique(grp_mat$o)) + 1):length(fulls)], dim(grp_mat)[1]/length(unique(grp_mat$o))/2))

## lists of all of the Y's
source(paste0(vrc01_code_dir, "makeDataAndFunctions.R")) # read in all of the data
ys_1 <- list(Y.sens.resis[!is.na(Y.sens.resis)], Y.cens[!is.na(Y.cens)])
ys_2 <- list(Y2.sens.resis[!is.na(Y2.sens.resis)], Y2.cens[!is.na(Y2.cens)])

## lists with copied Y's
ys_grp <- apply(grp_mat, 1, match_y, y1 = ys_1, y2 = ys_2, folds1 = folds[1:2], folds2 = folds[3:4], ord = c("sens.resis", "cens"))

## -------------------------------------
## variable importance analysis
## -------------------------------------
grp_indx <- min(which(grp_mat$s == 2)) - 1
start_time <- Sys.time()
for (i in 1:length(grp_nms)) {
    eval(parse(text = paste0("suppressWarnings(", grp_nms[[i]], " <- cv_vim(Y = ys_grp[[i]]$y, f1 = fulls_grp[[i]], f2 = red_grps[[i]], indx = rep(rep(1:num_grp, each = length(unique(full_mat$o))), 2)[i], V = 10, folds = ys_grp[[i]]$folds, type = args$vimp_measure, run_regression = FALSE, na.rm = TRUE, alpha = 1 - args$level, scale = 'identity'))")))
}
end_time <- Sys.time()
print(end_time - start_time)

## ---------------------------------------------------
## group importance; average, then combine
## ---------------------------------------------------
## average the importance across the two datasets
for (i in 1:length(nms_group)) {
  ## extract the feature set
  f_set_nms <- get_feature_list(grp_nms, i, grp_mat)
  ## average the results across datasets
  for (j in 1:length(unique(full_mat$o))) {
    ## logical vector for which to select
    logivec <- unlist(lapply(f_set_nms, function(x) grepl(unique(full_mat$o)[j], x)))
    ## comma separated vector of names
    nmvec <- paste(unlist(f_set_nms[logivec]), collapse = ", ")
    eval(parse(text = paste0("grp_", unique(full_mat$o)[j], "_set_", i,
                             " <- average_vim(", nmvec, ", weights = c(1/2, 1/2))")))
  }
}
## combine estimates for each outcome
for (i in 1:length(unique(grp_mat_chr$o))) {
    eval(parse(text = paste0("grp_", unique(grp_mat_chr$o)[i], "_avg",
                         " <- merge_vim(", paste(paste0("grp_", unique(grp_mat_chr$o)[i], "_set_", 1:num_grp), collapse = ", "),
                         ")")))
    eval(parse(text = paste0("grp_", unique(grp_mat_chr$o)[i], "_avg$s",
                         " <- unlist(lapply(strsplit(grp_", unique(grp_mat_chr$o)[i], "_avg$s, '_', fixed = TRUE), function(x) tail(x, n = 2)[1]))")))
    eval(parse(text = paste0("saveRDS(grp_", unique(grp_mat_chr$o)[i],   "_avg", ", file = paste0(output_dir, args$vimp_measure, '_vim_avg_grp_", unique(grp_mat_chr$o)[i], ".rds'))")))
}


## ----------------------------------------------------
## all estimates, for each outcome
## ----------------------------------------------------
## combine all of the estimates, for each outcome
for (j in 1:2) {
    for (i in 1:length(unique(grp_mat_chr$o))) {
      eval(parse(text = paste0("grp_", unique(grp_mat_chr$o)[i], "_", j,
                           " <- merge_vim(", paste(grp_nms[(1:grp_indx + grp_indx*(j-1))][grepl(unique(grp_mat_chr$o)[i], grp_nms[(1:grp_indx + grp_indx*(j-1))])], collapse = ", "),
                           ")")))
      eval(parse(text = paste0("saveRDS(grp_", unique(grp_mat_chr$o)[i], "_", j, ", file = paste0(output_dir, args$vimp_measure, '_vim_grp_", unique(grp_mat_chr$o)[i], "_", j, ".rds'))")))
    }
}
