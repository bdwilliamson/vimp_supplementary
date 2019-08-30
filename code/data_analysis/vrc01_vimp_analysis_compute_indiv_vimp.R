#!/usr/local/bin/Rscript
##########################################################################
##
## FILE: vrc01_vimp_analysis_compute_vimp.R
##
## PURPOSE: variable importance analysis of binary ic50 endpoints
##          from VRC01 CATNAP analysis, using deviance, accuracy, auc
##
##          this file computes variable importance estimates based on
##          estimated regression functions, and saves the main output
##          FOR INDIVIDUAL VARS ONLY
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
num_indiv <- 797

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

## load in reduced regression fitted values
ind_mat <- expand.grid(o = c("sens.resis", "cens"),
                     f = 1:num_indiv, s = 1:2, stringsAsFactors = FALSE)
ind_mat_chr <- ind_mat
ind_mat_chr$f <- as.character(ind_mat$f)
ind_mat_chr$s <- as.character(ind_mat$s)

## names of the fitted values
red_ind_list <- as.list(paste0(vrc01_results_dir, run, "/reduced_preds_",
                             ind_mat$o, "_dataset_", ind_mat$s, "_ind_",
                             ind_mat$f, "_", run,".rds"))


red_inds_rev <- lapply(red_ind_list, readRDS)
ind_nms <- apply(ind_mat_chr, 1, function(x) paste0("ind_", paste(x, collapse = "_")))

## reverse them all; the 10th-fold-fit is currently in the first position for each list
red_inds <- lapply(red_inds_rev, rev)

## create lists with the appropriate full fits, same length as both group and individual ones
num_reps <- dim(ind_mat)[1]/length(unique(ind_mat$o))/2

fulls_ind <- c(rep(fulls[1:length(unique(ind_mat$o))], num_reps),
                 rep(fulls[(length(unique(ind_mat$o)) + 1):length(fulls)], num_reps))

## lists of all of the Y's
source(paste0(vrc01_code_dir, "makeDataAndFunctions.R")) # read in all of the data
ys_1 <- list(Y.sens.resis[!is.na(Y.sens.resis)], Y.cens[!is.na(Y.cens)])
ys_2 <- list(Y2.sens.resis[!is.na(Y2.sens.resis)], Y2.cens[!is.na(Y2.cens)])

## lists with copied Y's
ys_ind <- apply(ind_mat, 1, match_y, y1 = ys_1, y2 = ys_2, folds1 = folds[1:2], folds2 = folds[3:4], ord = c("sens.resis", "cens"))

## names for individual
nms_ind <- unique(predictors)

## -------------------------------------
## variable importance analysis
## -------------------------------------
ind_indx <- min(which(ind_mat$s == 2)) - 1
start_time <- Sys.time()
for (i in 1:length(ind_nms)) {
    eval(parse(text = paste0("suppressWarnings(", ind_nms[[i]],
                               " <- cv_vim(Y = ys_ind[[i]]$y, f1 = fulls_ind[[i]], f2 = red_inds[[i]], indx = rep(rep(1:num_indiv, each = length(unique(full_mat$o))), 2)[i], V = 10, folds = ys_ind[[i]]$folds, type = args$vimp_measure, run_regression = FALSE, na.rm = TRUE, alpha = 1 - args$level, scale = 'identity'))")))
}
end_time <- Sys.time()
print(end_time - start_time)

## ---------------------------------------------------
## variable importance analysis: average, then combine
## ---------------------------------------------------
for (i in 1:length(nms_ind)) {
    ## extract the feature set
    f_set_nms <- get_feature_list(ind_nms, i, ind_mat)
    ## average the results across datasets
    for (j in 1:length(unique(full_mat$o))) {
      ## logical vector for which to select
      logivec <- unlist(lapply(f_set_nms, function(x) grepl(unique(full_mat$o)[j], x)))
      ## comma separated vector of names
      nmvec <- paste(unlist(f_set_nms[logivec]), collapse = ", ")
      if (i == 563) {
          eval(parse(text = paste0("ind_", unique(full_mat$o)[j], "_set_", i,
                               " <- average_vim(", nmvec, ", weights = c(1/2, 1/2))")))
        } else {
        eval(parse(text = paste0("ind_", unique(full_mat$o)[j], "_set_", i,
                             " <- average_vim(", nmvec, ", weights = c(1/2, 1/2))")))
      }
    }
}
## combine estimates for each outcome
for (i in 1:length(unique(ind_mat_chr$o))) {
    eval(parse(text = paste0("ind_", unique(ind_mat_chr$o)[i], "_avg",
                         " <- merge_vim(", paste(paste0("ind_", unique(ind_mat_chr$o)[i], "_set_", 1:num_indiv), collapse = ", "),
                         ")")))
    eval(parse(text = paste0("ind_", unique(ind_mat_chr$o)[i], "_avg$s",
                         " <- unlist(lapply(strsplit(ind_", unique(ind_mat_chr$o)[i], "_avg$s, '_', fixed = TRUE), function(x) tail(x, n = 2)[1]))")))
    eval(parse(text = paste0("saveRDS(ind_", unique(ind_mat_chr$o)[i],   "_avg", ", file = paste0(output_dir, args$vimp_measure, '_vim_avg_ind_", unique(ind_mat_chr$o)[i], ".rds'))")))
}

## ----------------------------------------------------
## all estimates, for each outcome
## ----------------------------------------------------
## combine all of the estimates, for each outcome
for (j in 1:2) {
    for (i in 1:length(unique(ind_mat_chr$o))) {
      vec <- ind_nms[(1:ind_indx + ind_indx*(j-1))][grepl(unique(ind_mat_chr$o)[i], ind_nms[(1:ind_indx + ind_indx*(j-1))])]

      eval(parse(text = paste0("ind_", unique(ind_mat_chr$o)[i], "_", j,
                           " <- merge_vim(", paste(vec, collapse = ", "),
                           ")")))
      eval(parse(text = paste0("saveRDS(ind_", unique(ind_mat_chr$o)[i], "_", j, ", file = paste0(output_dir, args$vimp_measure, '_vim_ind_", unique(ind_mat_chr$o)[i], "_", j, ".rds'))")))
    }
}
