#!/usr/local/bin/Rscript

# get variable importance estimates

# ---------------------------------------------------------------------------
# Set up args, variables, functions
# ---------------------------------------------------------------------------

# load libraries
library("SuperLearner")
library("dplyr")
library("argparse")
library("here")

lib_paths <- .libPaths()
# require vimp@v2.2.3 or higher
if (is.na(Sys.getenv("RSTUDIO", unset = NA))) {
    library("vimp", lib.loc = lib_paths[2])
    code_dir <- ""
} else {
    library("vimp")
    code_dir <- "code/data_analysis"
}

source(here(code_dir, "utils.R"))
source(here(code_dir, "get_variable_groups.R"))
source(here(code_dir, "super_learner_libraries.R"))

parser <- ArgumentParser()
parser$add_argument("--outcome", default = c("sens50", "sens80", "cens"), 
                    help = paste0("outcome string: ",
                                  "cens for right-censored IC-50,",
                                  "sens50 for IC-50 < 1,", 
                                  "sens80 for IC-80 < 1.")
)
parser$add_argument("--ind-type", default = "sitewise",
                    help = "should we do site-wise or residue-wise individual importance?")
parser$add_argument("--reduce-groups", action = "store_true", help = "should we run only for CD4 binding sites, or for all groups?")
parser$add_argument("--reduce-covs", action = "store_true", help = "should we run only for only a subset of the covariates?")
parser$add_argument("--reduce-library", action = "store_true", help = "should we run with a small learner library, or with all candidate learners?")
parser$add_argument("--vimp-measure", default = "auc", help = "type of variable importance measure to use")
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

# ---------------------------------------------------------------------------
# get variable importance!
# ---------------------------------------------------------------------------
set.seed(4747)
for (i in seq_len(length(args$outcome))) {
    this_outcome_name <- get_outcome(args$outcome[i])
    this_analysis_dataset <- initial_analysis_dataset[, c(this_outcome_name, pred_names)]
    this_analysis_dataset_cc <- this_analysis_dataset[complete.cases(this_analysis_dataset), ]
    y <- this_analysis_dataset_cc[, this_outcome_name]
    # create output list
    eval(parse(
        text = paste0(this_outcome_name,
                      '_cv_vimp_lst <- make_vimp_list(var_groups, all_var_inds)')
    ))
    # load full fit corresponding to this outcome
    full_cv_fit <- readRDS(here(code_dir, paste0("slfits/cvfit_", this_outcome_name, ".rds")))
    cross_fitting_folds <- vimp::get_cv_sl_folds(full_cv_fit$folds)
    ss_V <- length(unique(cross_fitting_folds))
    V <- ss_V / 2
    sample_splitting_folds <- vimp::make_folds(
        y = unique(cross_fitting_folds), V = 2
    )
    full_cv_preds <- vimp::extract_sampled_split_predictions(
        cvsl_obj = full_cv_fit, sample_splitting = TRUE,
        sample_splitting_folds = sample_splitting_folds, full = TRUE
    )
    # full_fit <- readRDS(here(code_dir, paste0("slfits/fitted_", this_outcome_name, ".rds")))
    # load geog-only fit corresponding to this outcome
    geog_cv_fit <- readRDS(here(code_dir, paste0("slfits/cvfit_", this_outcome_name, "_geog.rds")))
    geog_cv_preds <- vimp::extract_sampled_split_predictions(
        cvsl_obj = geog_cv_fit, sample_splitting = TRUE,
        sample_splitting_folds = sample_splitting_folds, full = FALSE
    )
    # geog_fit <- readRDS(here(code_dir, paste0("slfits/fitted_", this_outcome_name, "_geog.rds")))
    # group variable importance
    for (j in 1:length(var_groups)) {
        this_group_name <- names(var_groups)[j]
        # cv
        cond_cv_fit <- readRDS(here(code_dir, paste0("slfits/cvfit_",
                                                     this_outcome_name, "_cond_",
                                                     this_group_name, ".rds")))
        marg_cv_fit <- readRDS(here(code_dir, paste0("slfits/cvfit_",
                                                     this_outcome_name, "_marg_",
                                                     this_group_name, ".rds")))
        cond_cv_preds <- vimp::extract_sampled_split_predictions(
            cvsl_obj = cond_cv_fit, sample_splitting = TRUE,
            sample_splitting_folds = sample_splitting_folds, full = FALSE
        )
        marg_cv_preds <- vimp::extract_sampled_split_predictions(
            cvsl_obj = marg_cv_fit, sample_splitting = TRUE,
            sample_splitting_folds = sample_splitting_folds, full = TRUE
        )
        # cond_fit <- readRDS(here(code_dir, paste0("slfits/fitted_",
        #                                           this_outcome_name, "_cond_",
        #                                           this_group_name, ".rds")))
        # marg_fit <- readRDS(here(code_dir, paste0("slfits/fitted_",
        #                                           this_outcome_name, "_marg_",
        #                                           this_group_name, ".rds")))
        # get conditional, cv vimp
        suppressWarnings(eval(parse(text = paste0(
            this_outcome_name, "_vim_cond_", this_group_name,
            " <- cv_vim(Y = y, cross_fitted_f1 = full_cv_preds, cross_fitted_f2 = cond_cv_preds,
             indx = which(pred_names %in% var_groups[[j]]),
             delta = 0, V = V, run_regression = FALSE, sample_splitting = TRUE,
             cross_fitting_folds = cross_fitting_folds, sample_splitting_folds = sample_splitting_folds,
             type = args$vimp_measure)"
        ))))
        # get marginal, cv vimp
        suppressWarnings(eval(parse(text = paste0(
            this_outcome_name, "_vim_marg_", this_group_name,
            " <- cv_vim(Y = y, cross_fitted_f1 = marg_cv_preds, cross_fitted_f2 = geog_cv_preds,
             indx = which(pred_names %in% var_groups[[j]]),
             delta = 0, V = 5, run_regression = FALSE, sample_splitting = TRUE,
             cross_fitting_folds = cross_fitting_folds, sample_splitting_folds = sample_splitting_folds,
             type = args$vimp_measure)"
        ))))
    }
    # merge together
    eval(parse(text = paste0(this_outcome_name,
                             "_cv_vimp_lst$conditional <- merge_vim(",
                             paste(paste0(this_outcome_name, "_vim_cond_",
                                          names(var_groups)), collapse = ", "), ")")))
    eval(parse(text = paste0(this_outcome_name,
                             "_cv_vimp_lst$marginal <- merge_vim(",
                             paste(paste0(this_outcome_name, "_vim_marg_",
                                          names(var_groups)[-length(names(var_groups))]), 
                                   collapse = ", "), ")")))
    # individual variable importance
    for (j in 1:length(all_var_inds)) {
        this_var_name <- names(all_var_inds)[j]
        # cv
        ind_cv_fit <- readRDS(here(code_dir, paste0("slfits/cvfit_", 
                                                    this_outcome_name, "_marg_",
                                                    this_var_name, ".rds")))
        ind_cv_preds <- vimp::extract_sampled_split_predictions(
            cvsl_obj = ind_cv_fit, sample_splitting = TRUE,
            sample_splitting_folds = sample_splitting_folds, full = TRUE
        )
        # ind_fit <- readRDS(here(code_dir, paste0("slfits/fitted_", 
        #                                          this_outcome_name, "_marg_",
        #                                          this_var_name, ".rds")))
        # get individual, cv vimp
        suppressWarnings(eval(parse(text = paste0(
            this_outcome_name, "_vim_marg_", this_var_name,
            " <- cv_vim(Y = y, cross_fitted_f1 = ind_cv_preds, cross_fitted_f2 = geog_cv_preds,
             indx = which(pred_names %in% all_var_inds[[j]]),
             delta = 0, V = 5, run_regression = FALSE, sample_splitting = TRUE,
             cross_fitting_folds = cross_fitting_folds, sample_splitting_folds = sample_splitting_folds,
             type = args$vimp_measure)"
        ))))
    }
    # merge together
    eval(parse(text = paste0(this_outcome_name,
                             "_cv_vimp_lst$individual <- merge_vim(",
                             paste(paste0(this_outcome_name, "_vim_marg_",
                                          names(all_var_inds)), collapse = ", "), ")")))
    # save them off
    eval(parse(text = paste0(
        "saveRDS(", this_outcome_name, "_cv_vimp_lst, file = here(code_dir, 'slfits/",
        paste0(this_outcome_name, "_", args$vimp_measure, "_cv_vimp"), ".rds'))"
    )))
}
