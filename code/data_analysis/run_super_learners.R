#!/usr/local/bin/Rscript

## run super learners for VRC01 CATNAP analysis

## ----------------------------------------------------------------
## set up
## ----------------------------------------------------------------
library("SuperLearner")
library("argparse")
source("../code/get_variable_groups.R")
source("../code/super_learner_libraries.R")

parser <- ArgumentParser()
parser$add_argument("--outcome", default = "cens", help = "outcome string")
parser$add_argument("-redgrps", "--reduce-groups", action = "store_true", help = "should we run only for CD4 binding sites, or for all groups?")
parser$add_argument("-redcovs", "--reduce-covs", action = "store_true", help = "should we run only for only a subset of the covariates?")
parser$add_argument("-redlib", "--reduce-library", action = "store_true", help = "should we run with a small learner library, or with all candidate learners?")
parser$add_argument("-rungeog", "--run-geog", action = "store_false", help = "should we run geog sls?")
parser$add_argument("-runindi", "--run-indiv", action = "store_false", help = "should we run individual sls?")
parser$add_argument("-runfull", "--run-full", action = "store_false", help = "should we run the full sl?")
parser$add_argument("-rungrps", "--run-grps", action = "store_false", help = "should we run the minus group sls?")
args <- parser$parse_args()
print(args)

## read in the data
dat <- readRDS("../code/data/analysis_data.rds")

## if reduce_library, run with small lib
if (args$reduce_library) {
    learner_lib <- library_reduced
} else {
    learner_lib <- library_full
}

## get names of predictors
all_var_groups <- get_variable_groups(dat)
## get variable groups
all_geog_vars <- all_var_groups$geog
## get all predictor names
pred_names <- unique(do.call(c, all_var_groups))
# if reduce_groups, only run on the cd4 binding site
if (args$reduce_groups) {
    var_groups <- all_var_groups[1]
} else {
    var_groups <- all_var_groups
}
# get names of outcomes
all_outcome_names <- c("ic50.censored", "binding.dichotomous.sens.resis")
outcome_names <- all_outcome_names[grepl(args$outcome, all_outcome_names)]
dat$ic50.censored <- as.numeric(dat$ic50.censored)
## specify the individual variables
if (args$reduce_covs) {
    num_covs <- 10
    var_inds <- pred_names[!grepl("geog", pred_names)][1:(num_covs - length(all_geog_vars))]
} else {
    num_covs <- length(pred_names) - length(all_geog_vars)
    var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]
}
V <- 10
# determine SL options based on outcome name
get_sl_options <- function(outcome_name) {
    if (grepl("cens", outcome_name) | grepl("sens.resis", outcome_name)) {
        sl_fam <- "binomial"
        cv_ctrl_lst <- list(V = V, stratifyCV = TRUE)
        sl_method <- "tmp_method.CC_nloglik"
    } else {
        sl_fam <- "gaussian"
        cv_ctrl_lst <- list(V = V)
        sl_method <- "tmp_method.CC_LS"
    }
    return(list(fam = sl_fam, ctrl = cv_ctrl_lst, method = sl_method))
}


make_folds <- function(y, V, stratified = TRUE) {
  if (stratified) {
    y_1 <- y == 1
    y_0 <- y == 0
    folds_1 <- rep(seq_len(V), length = sum(y_1))
    folds_1 <- sample(folds_1)
    folds_0 <- rep(seq_len(V), length = sum(y_0))
    folds_0 <- sample(folds_0)
    folds <- vector("numeric", length(draw$y_cat))
    folds[y_1] <- folds_1
    folds[y_0] <- folds_0
  } else {
    folds <- rep(seq_len(V), length = length(y))
    folds <- sample(folds)
  }
  return(folds)
}

#' function to run cv super learner on a single outcome
#' @param outcome_name String name of outcome
#' @param pred_names Vector of string names of predictor variables
#' @param save_dir name of directory to save results to
#' @param cv_fit_name name of CV fits (defaults to cvfit_<outcome_name>.rds)
#' @param full_fit is it the full fit (TRUE) or a reduced fit (FALSE)?
#' @param save_full_object Flag for whether or not to save the full fitted object, or just the fitted values
#' @param ... other arguments to super learner
sl_one_outcome <- function(outcome_name, pred_names, save_dir = "../code/slfits/",  cv_fit_name = paste0("cvfit_", outcome_name, ".rds"), full_fit = TRUE, save_full_object = TRUE, reduce_covs = FALSE, ...){

    if (full_fit) {
        outcome_bool <- dat$dataset == 1
    } else {
        outcome_bool <- dat$dataset == 2
    }
    newdat <- subset(dat, outcome_bool)

    pred <- newdat[ , pred_names]

    if (reduce_covs) {
        pred <- pred[, 1:num_covs]
    }

    y <- newdat[, outcome_name]
    Y <- y[!is.na(y)]
    preds <- pred[!is.na(y), ]


    cv_fit <- CV.SuperLearner(Y = Y, X = preds, ...)
    if (save_full_object) {
        saveRDS(cv_fit, file = paste0(save_dir, cv_fit_name))
    }
    saveRDS(cv_fit$SL.predict, file = paste0(save_dir, gsub(".RData", ".rds", gsub("cvfit_", "cvfitted_", cv_fit_name))))
    saveRDS(cv_fit$folds, file = paste0(save_dir, gsub("cvfitted_", "cvfolds_", gsub(".RData", ".rds", gsub("cvfit_", "cvfitted_", cv_fit_name)))))

    return(invisible(NULL))
}
## create directory if necessary
if (!dir.exists("../code/slfits/")) {
    dir.create("../code/slfits/")
}
## ----------------------------------------------------------------------------
## (1) run full super learners for each outcome (unless reduce_outcomes = TRUE)
## ----------------------------------------------------------------------------
set.seed(474747)
if (args$run_full) {
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        cat("Running full SL: ", this_outcome_name, "\n")
        sl_opts <- get_sl_options(this_outcome_name)
        sl_fit_i <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names, reduce_covs = args$reduce_covs, family = sl_opts$fam, SL.library = learner_lib, cvControl = sl_opts$ctrl, method = sl_opts$method)
    }
}


## ----------------------------------------------------------------------------
## (2)+(3) run super learners for each outcome (unless reduce_outcomes = TRUE)
##         on (2) reduced set of features defined by removing group of interest
##         and (3) group of interest + geographic confounders
## ----------------------------------------------------------------------------
## run super learners on pre-defined groups
if (args$run_grps) {
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        cat("Running reduced SL: ", this_outcome_name, "\n")
        sl_opts <- get_sl_options(this_outcome_name)
        ## read in the full folds, for this group's cv folds
        # full_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, ".rds"))
        # sl_opts$ctrl$validRows <- full_folds
        for (j in 1:length(var_groups)) {
            if (length(var_groups[j]) != 0) {
                this_group_name <- names(var_groups)[j]
                cat("Group: ", this_group_name, "\n")
                ## fit based on removing group of interest
                sl_fit_ij <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[!(pred_names %in% var_groups[[j]])], cv_fit_name = paste0("cvfitted_", this_outcome_name, "_minus_", this_group_name, ".rds"), reduce_covs = args$reduce_covs, family = sl_opts$fam, SL.library = learner_lib, cvControl = sl_opts$ctrl, method = sl_opts$method, save_full_object = FALSE, full_fit = FALSE)
                ## fit based on only group of interest + geographic confounders
                sl_fit_marginal_ij <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% var_groups[[j]]) | (pred_names %in% all_geog_vars)], cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"), reduce_covs = args$reduce_covs, family = sl_opts$fam, SL.library = learner_lib, cvControl = sl_opts$ctrl, method = sl_opts$method, save_full_object = FALSE, full_fit = FALSE)
            }
        }
    }
}
## ----------------------------------------------------------------------------
## (4) run super learners for each outcome (unless reduce_outcomes = TRUE)
##     on each individual feature + confounders
## ----------------------------------------------------------------------------
if (args$run_indiv) {
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        cat("Running individual SL: ", this_outcome_name, "\n")
        sl_opts <- get_sl_options(this_outcome_name)
        # full_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, ".rds"))
        # sl_opts$ctrl$validRows <- full_folds
        for (j in 1:length(var_inds)) {
            this_var_name <- var_inds[j]
            ## fit SL of only this variable plus geographic confounders
            sl_fit_ij <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% var_inds[j]) | (pred_names %in% all_geog_vars)], cv_fit_name = paste0("cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"), reduce_covs = FALSE, save_full_object = FALSE, family = sl_opts$fam, SL.library = learner_lib, cvControl = sl_opts$ctrl, method = sl_opts$method, full_fit = FALSE)
        }
    }
}

## ----------------------------------------------------------------------------
## (5) run super learners for each outcome (unless reduce_outcomes = TRUE)
##     on only confounders
## ----------------------------------------------------------------------------
if (args$run_geog) {
    for (i in 1:length(outcome_names)) {
        this_outcome_name <- outcome_names[i]
        cat("Running geog SL: ", this_outcome_name, "\n")
        sl_opts <- get_sl_options(this_outcome_name)
        # full_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, ".rds"))
        # sl_opts$ctrl$validRows <- full_folds
        sl_geog_i <- sl_one_outcome(outcome_name = this_outcome_name, pred_names = pred_names[(pred_names %in% all_geog_vars)], cv_fit_name = paste0("cvfitted_", this_outcome_name, "_geog.rds"), reduce_covs = FALSE, save_full_object = FALSE, family = sl_opts$fam, SL.library = learner_lib, cvControl = sl_opts$ctrl, method = sl_opts$method, full_fit = TRUE)
    }
}
