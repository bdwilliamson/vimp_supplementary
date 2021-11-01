#!/usr/local/bin/Rscript

# run super learners for VRC01 CATNAP analysis

# ----------------------------------------------------------------
# set up
# ----------------------------------------------------------------
library("SuperLearner")
library("argparse")
library("here")

if (is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    code_dir <- ""
} else {
    code_dir <- "code"
}
source(here(code_dir, "get_variable_groups.R"))
source(here(code_dir, "super_learner_libraries.R"))
source(here(code_dir, "utils.R"))

parser <- ArgumentParser()
parser$add_argument("--outcome", default = "cens", 
                    help = paste0("outcome string: ",
                                  "cens for right-censored IC-50,",
                                  "sens50 for IC-50 < 1,", 
                                  "sens80 for IC-80 < 1.")
                    )
parser$add_argument("--ind-type", default = "sitewise",
                    help = "should we do site-wise or residue-wise individual importance?")
# the next args **only store the action if they are present**;
# the default is the opposite of the action
parser$add_argument("-redgrps", "--reduce-groups", action = "store_true",
                    help = "should we run only for CD4 binding sites, or for all groups?")
parser$add_argument("-redcovs", "--reduce-covs", action = "store_true",
                    help = "should we run only for only a subset of the covariates?")
parser$add_argument("-redlib", "--reduce-library", action = "store_true",
                    help = "should we run with a small learner library, or with all candidate learners?")
parser$add_argument("-dontrungeog", "--dontrun-geog", action = "store_true",
                    help = "should we not run geog sls?")
parser$add_argument("-dontrunindi", "--dontrun-indiv", action = "store_true",
                    help = "should we not run individual sls?")
parser$add_argument("-dontrunfull", "--dontrun-full", action = "store_true",
                    help = "should we not run the full sl?")
parser$add_argument("-dontrungrps", "--dontrun-grps", action = "store_true",
                    help = "should we not run the minus group sls?")
args <- parser$parse_args()
print(args)

# read in the data
dat <- readRDS(here(code_dir, "data/analysis_data.rds"))

# if reduce_library, run with small lib
if (args$reduce_library) {
    learner_lib <- library_reduced
} else {
    learner_lib <- library_full
}
print(learner_lib)

# specify the individual variables and groups
# get names of predictors
all_var_groups <- get_variable_groups(dat)
# get variable groups
all_geog_vars <- all_var_groups$geog
# get all predictor names
pred_names <- unique(do.call(c, all_var_groups))
if (args$reduce_covs) {
    num_covs <- 10
} else {
    num_covs <- length(pred_names)
}
# get all individual predictors
all_var_inds <- get_individual_features(pred_names[!grepl("geog", pred_names)][1:num_covs],
                                        args$ind_type)

# if reduce_groups, only run on the cd4 binding site
if (args$reduce_groups) {
    var_groups <- all_var_groups[1]
} else {
    var_groups <- all_var_groups
}
# the only possible outcomes for these analyses are: IC50 censored ("ic50.censored"),
# sensitive/resistant ("binding.dichotomous.sens.resis"), IC50 < 1, and IC80 < 1
initial_analysis_dataset <- dat[, c("ic50.censored", "binding.dichotomous.sens.resis",
                                    "sens50", "sens80", pred_names)]
# specify the number of folds; note that this allows 5-fold cross-fitted VIM estimation
V <- 10

# create directory if necessary
if (!dir.exists(here(code_dir, "slfits"))) {
    dir.create(here(code_dir, "slfits"))
}
# ----------------------------------------------------------------------------
# (1) run full super learners for each outcome (unless reduce_outcomes = TRUE)
# ----------------------------------------------------------------------------
set.seed(474747)
if (!args$dontrun_full) {
    for (i in seq_len(length(args$outcome))) {
        this_outcome_name <- get_outcome(args$outcome[i])
        this_analysis_dataset <- initial_analysis_dataset[, c(this_outcome_name, pred_names)]
        this_analysis_dataset_cc <- this_analysis_dataset[complete.cases(this_analysis_dataset), ]
        cat("Running full SL: ", this_outcome_name, "\n")
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        sl_fit_i <- sl_one_outcome(
          dat = this_analysis_dataset_cc, save_dir = here(code_dir, "slfits/"),
          outcome_name = this_outcome_name, pred_names = pred_names,
          reduce_covs = args$reduce_covs, family = sl_opts$fam, SL.library = learner_lib,
          cvControl = sl_opts$ctrl, method = sl_opts$method, innerCvControl = list(list(V = 5))
        )
    }
}

# ----------------------------------------------------------------------------
# (2)+(3) run super learners for each outcome (unless reduce_outcomes = TRUE)
#         on (2) reduced set of features defined by removing group of interest
#         and (3) group of interest + geographic confounders
# ----------------------------------------------------------------------------
# run super learners on pre-defined groups
if (!args$dontrun_grps) {
    for (i in seq_len(length(args$outcome))) {
        this_outcome_name <- get_outcome(args$outcome[i])
        this_analysis_dataset <- initial_analysis_dataset[, c(this_outcome_name, pred_names)]
        this_analysis_dataset_cc <- this_analysis_dataset[complete.cases(this_analysis_dataset), ]
        cat("Running reduced SL: ", this_outcome_name, "\n")
        cross_fitting_folds <- readRDS(paste0("slfits/cvfolds_", this_outcome_name, ".rds"))
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        sl_opts$ctrl$validRows <- cross_fitting_folds
        for (j in 1:length(var_groups)) {
            if (length(var_groups[j]) != 0) {
                this_group_name <- names(var_groups)[j]
                cat("Group: ", this_group_name, "\n")
                # fit based on removing group of interest
                sl_fit_cond_ij <- sl_one_outcome(
                    dat = this_analysis_dataset_cc, save_dir = here(code_dir, "slfits/"),
                    outcome_name = this_outcome_name,
                    pred_names = pred_names[!(pred_names %in% var_groups[[j]])],
                    cv_fit_name = paste0("cvfit_", this_outcome_name, "_cond_", this_group_name, ".rds"),
                    reduce_covs = args$reduce_covs, family = sl_opts$fam,
                    SL.library = learner_lib, cvControl = sl_opts$ctrl,
                    method = sl_opts$method, save_full_object = TRUE,
                    innerCvControl = list(list(V = 5))
                )
                # fit based on only group of interest + geographic confounders
                sl_fit_marg_ij <- sl_one_outcome(
                    dat = this_analysis_dataset_cc, save_dir = here(code_dir, "slfits/"),
                    outcome_name = this_outcome_name,
                    pred_names = pred_names[(pred_names %in% var_groups[[j]]) | (pred_names %in% all_geog_vars)],
                    cv_fit_name = paste0("cvfit_", this_outcome_name, "_marg_", this_group_name, ".rds"),
                    reduce_covs = args$reduce_covs, family = sl_opts$fam,
                    SL.library = learner_lib, cvControl = sl_opts$ctrl,
                    method = sl_opts$method, save_full_object = TRUE,
                    innerCvControl = list(list(V = 5))
                )
            }
        }
    }
}
# ----------------------------------------------------------------------------
# (4) run super learners for each outcome (unless reduce_outcomes = TRUE)
#     on each individual feature + confounders
# ----------------------------------------------------------------------------
if (!args$dontrun_indiv) {
    for (i in seq_len(length(args$outcome))) {
        this_outcome_name <- get_outcome(args$outcome[i])
        this_analysis_dataset <- initial_analysis_dataset[, c(this_outcome_name, pred_names)]
        this_analysis_dataset_cc <- this_analysis_dataset[complete.cases(this_analysis_dataset), ]
        cat("Running individual SL: ", this_outcome_name, "\n")
        cross_fitting_folds <- readRDS(paste0("slfits/cvfolds_", this_outcome_name, ".rds"))
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        sl_opts$ctrl$validRows <- cross_fitting_folds
        for (j in 1:length(all_var_inds)) {
            this_var_name <- names(all_var_inds)[j]
            # fit SL of only this variable plus geographic confounders
            sl_fit_ind_ij <- sl_one_outcome(
                dat = this_analysis_dataset_cc, save_dir = here(code_dir, "slfits/"),
                outcome_name = this_outcome_name,
                pred_names = pred_names[(pred_names %in% all_var_inds[[j]]) | (pred_names %in% all_geog_vars)],
                cv_fit_name = paste0("cvfit_", this_outcome_name, "_marg_", this_var_name, ".rds"),
                reduce_covs = FALSE, save_full_object = TRUE, family = sl_opts$fam,
                SL.library = learner_lib, cvControl = sl_opts$ctrl,
                method = sl_opts$method, innerCvControl = list(list(V = 5))
            )
        }
    }
}

# ----------------------------------------------------------------------------
# (5) run super learners for each outcome (unless reduce_outcomes = TRUE)
#     on only confounders
# ----------------------------------------------------------------------------
if (!args$dontrun_geog) {
    for (i in seq_len(length(args$outcome))) {
        this_outcome_name <- get_outcome(args$outcome[i])
        this_analysis_dataset <- initial_analysis_dataset[, c(this_outcome_name, pred_names)]
        this_analysis_dataset_cc <- this_analysis_dataset[complete.cases(this_analysis_dataset), ]
        cat("Running geog SL: ", this_outcome_name, "\n")
        cross_fitting_folds <- readRDS(paste0("slfits/cvfolds_", this_outcome_name, ".rds"))
        sl_opts <- get_sl_options(this_outcome_name, V = V)
        sl_opts$ctrl$validRows <- cross_fitting_folds
        sl_geog_i <- sl_one_outcome(
            dat = this_analysis_dataset_cc, save_dir = here(code_dir, "slfits/"),
            outcome_name = this_outcome_name,
            pred_names = pred_names[(pred_names %in% all_geog_vars)],
            cv_fit_name = paste0("cvfit_", this_outcome_name, "_geog.rds"),
            reduce_covs = FALSE, save_full_object = TRUE, family = sl_opts$fam,
            SL.library = learner_lib, cvControl = sl_opts$ctrl,
            method = sl_opts$method, innerCvControl = list(list(V = 5))
        )
    }
}
