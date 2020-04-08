#!/usr/local/bin/Rscript

## get variable importance estimates

## ---------------------------------------------------------------------------
## Set up args, variables, functions
## ---------------------------------------------------------------------------

# load libraries
library("SuperLearner")
library("vimp")
library("dplyr")
library("argparse")
source("../code/get_variable_groups.R")
source("../code/super_learner_libraries.R")

parser <- ArgumentParser()
parser$add_argument("--outcome", default = "cens", help = "outcome string")
parser$add_argument("--reduce-groups", action = "store_true", help = "should we run only for CD4 binding sites, or for all groups?")
parser$add_argument("--reduce-covs", action = "store_true", help = "should we run only for only a subset of the covariates?")
parser$add_argument("--reduce-library", action = "store_true", help = "should we run with a small learner library, or with all candidate learners?")
parser$add_argument("--vimp-measure", default = "auc", help = "type of variable importance measure to use")
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
## get outer folds
outer_folds <- dat$dataset
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

make_vimp_list <- function(var_groups, var_inds) {
    list_names <- c("conditional", "marginal", "individual")
    lst <- sapply(list_names, function(x) NULL, simplify = FALSE)
    return(lst)
}
get_cv_folds <- function(folds_lst) {
    V <- length(folds_lst)
    v_lst <- sapply(1:V, function(s) rep(s, length(folds_lst[[s]])), simplify = FALSE)
    joint_lst <- mapply(list, v_lst, folds_lst, SIMPLIFY = FALSE)
    folds_mat <- do.call(rbind, lapply(joint_lst, function(x) cbind(x[[1]], x[[2]])))
    folds <- folds_mat[order(folds_mat[, 2]), 1]
    return(folds)
}
make_cv_lists <- function(folds_lst, full_vec, redu_vec) {
    folds <- get_cv_folds(folds_lst)
    ## make lists of the fitted values
    full_lst <- lapply(as.list(1:length(unique(folds))), function(x) full_vec[folds == x])
    redu_lst <- lapply(as.list(1:length(unique(folds))), function(x) redu_vec[folds == x])
    return(list(folds = folds, full_lst = full_lst, redu_lst = redu_lst))
}


## ---------------------------------------------------------------------------
## get variable importance!
## ---------------------------------------------------------------------------
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    ## create output list
    eval(parse(text = paste0(this_outcome_name, '_cv_vimp_lst <- make_vimp_list(var_groups, var_inds)')))
    ## load full fit corresponding to this outcome
    full_cv_fit <- readRDS(paste0("../code/slfits/cvfitted_", this_outcome_name, ".rds"))
    full_cv_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, ".rds"))
    full_cv_folds_vec <- get_cv_folds(full_cv_folds)
    full_cv_fit_lst <- lapply(as.list(1:length(unique(full_cv_folds_vec))), function(x) full_cv_fit[full_cv_folds_vec == x])
    ## load geog-only fit corresponding to this outcome
    geog_cv_fit <- readRDS(paste0("../code/slfits/cvfitted_", this_outcome_name, "_geog.rds"))
    geog_cv_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, "_geog.rds"))
    geog_cv_folds_vec <- get_cv_folds(geog_cv_folds)
    geog_cv_fit_lst <- lapply(as.list(1:length(unique(geog_cv_folds_vec))), function(x) geog_cv_fit[geog_cv_folds_vec == x])
    ## group variable importance
    for (j in 1:length(var_groups)) {
        this_group_name <- names(var_groups)[j]
        ## cv
        cond_cv_fit <- readRDS(paste0("../code/slfits/cvfitted_", this_outcome_name, "_minus_", this_group_name, ".rds"))
        marg_cv_fit <- readRDS(paste0("../code/slfits/cvfitted_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
        cond_cv_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, "_minus_", this_group_name, ".rds"))
        marg_cv_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, "_marginal_", this_group_name, ".rds"))
        cond_cv_folds_vec <- get_cv_folds(cond_cv_folds)
        marg_cv_folds_vec <- get_cv_folds(marg_cv_folds)
        cond_cv_fit_lst <- lapply(as.list(1:length(unique(cond_cv_folds_vec))), function(x) cond_cv_fit[cond_cv_folds_vec == x])
        marg_cv_fit_lst <- lapply(as.list(1:length(unique(marg_cv_folds_vec))), function(x) marg_cv_fit[marg_cv_folds_vec == x])
        ## get conditional, cv vimp
        # set.seed(474747)
        cond_folds <- list(outer_folds = outer_folds, inner_folds = list(inner_folds_1 = full_cv_folds_vec, inner_folds_2 = cond_cv_folds_vec))
        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_cond_", this_group_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = full_cv_fit_lst, f2 = cond_cv_fit_lst, indx = which(pred_names %in% var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = args$vimp_measure, folds = cond_folds, V = 10, na.rm = TRUE, scale = 'identity')"))))
        ## get marginal, cv vimp
        # set.seed(474747)
        marg_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2, inner_folds = list(inner_folds_1 = marg_cv_folds_vec, inner_folds_2 = geog_cv_folds_vec))
        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_marg_", this_group_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = marg_cv_fit_lst, f2 = geog_cv_fit_lst, indx = which(pred_names %in% var_groups[[j]]), run_regression = FALSE, alpha = 0.05, delta = 0, type = args$vimp_measure, folds = marg_folds, V = 10, na.rm = TRUE, scale = 'identity')"))))
    }
    ## merge together
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$conditional <- merge_vim(", paste(paste0(this_outcome_name, "_cv_cond_", names(var_groups)), collapse = ", "), ")")))
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$marginal <- merge_vim(", paste(paste0(this_outcome_name, "_cv_marg_", names(var_groups)), collapse = ", "), ")")))
    ## individual variable importance
    for (j in 1:length(var_inds)) {
        this_var_name <- var_inds[j]
        ## cv
        indi_cv_fit <- readRDS(paste0("../code/slfits/cvfitted_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
        indi_cv_folds <- readRDS(paste0("../code/slfits/cvfolds_", this_outcome_name, "_marginal_", this_var_name, ".rds"))
        indi_cv_folds_vec <- get_cv_folds(indi_cv_folds)
        indi_cv_fit_lst <- lapply(as.list(1:length(unique(indi_cv_folds_vec))), function(x) indi_cv_fit[indi_cv_folds_vec == x])
        ## get individual, cv vimp
        # set.seed(474747)
        indi_folds <- list(outer_folds = (-1) * (outer_folds - 1) + 2, inner_folds = list(inner_folds_1 = indi_cv_folds_vec, inner_folds_2 = geog_cv_folds_vec))
        suppressWarnings(eval(parse(text = paste0(this_outcome_name, "_cv_marg_", this_var_name, " <- vimp::cv_vim(Y = dat[, this_outcome_name], f1 = indi_cv_fit_lst, f2 = geog_cv_fit_lst, indx = which(pred_names %in% this_var_name), run_regression = FALSE, alpha = 0.05, delta = 0, type = args$vimp_measure, folds = indi_folds, V = 10, na.rm = TRUE, scale = 'identity')"))))
    }
    ## merge together
    eval(parse(text = paste0(this_outcome_name, "_cv_vimp_lst$individual <- merge_vim(", paste(paste0(this_outcome_name, "_cv_marg_", var_inds), collapse = ", "), ")")))
    ## save them off
    eval(parse(text = paste0("saveRDS(", this_outcome_name, "_cv_vimp_lst, file = '../code/slfits/", paste0(this_outcome_name, "_", args$vimp_measure, "_cv_vimp"), ".rds')")))
}
