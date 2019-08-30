#!/usr/local/bin/Rscript
## ------------------------------------------------
## FILE: sim_binary_bivariate.R
## CREATED: 06 March 2019 by Brian Williamson
## PURPOSE: general-purpose simulation setup for
##          binary outcome, bivariate predictor
## ------------------------------------------------

## ---------------------------------------------
## load required libraries and functions
## ---------------------------------------------
library("SuperLearner")
# only run this if necessary to update package
# devtools::install_github("bdwilliamson/vimp@v2.0.0")
library("vimp")
library("methods")
library("argparse")
library("xgboost")
library("gam")
library("dplyr")

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  source("lowdim/R/sim_binary_bivariate_data.R")
  source("lowdim/R/sim_binary_bivariate_ests.R")
  source("lowdim/R/run_sim_binary_bivariate_once.R")
} else {
  source("sim_binary_bivariate_data.R")
  source("sim_binary_bivariate_ests.R")
  source("run_sim_binary_bivariate_once.R")
}

## ---------------------------------------------
## pull in command-line arguments,
## set up the simulation
## ---------------------------------------------
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "bivariate_loss",
                    help = "the name of the simulation")
parser$add_argument("--risk-type", default = "expected_loss",
                    help = "whether or not the risk is an expected loss")
parser$add_argument("--vimp-measure", nargs = "+", default = "deviance",
                    help = "the measure of variable importance")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 50,
                    help = "number of replicates for each job")
parser$add_argument("--b", type = "double", default = 1000,
                    help = "number of bootstrap replicates")
parser$add_argument("--cv", type = "double", default = 1,
                    help = "run cross-validated estimators?")
args <- parser$parse_args()
print(args)

## true parameter values
p <- 0.6
if (grepl("null", args$sim_name)) {
    truths <- readRDS("truths_binary_bivariate_null.rds")
    mu_0 <- matrix(0, ncol = 1)
    mu_1 <- matrix(1.5, ncol = 1)
    Sigma <- diag(1, nrow = 1)
} else {
    truths <- readRDS("truths_binary_bivariate.rds")
    mu_0 <- matrix(c(0, 0), nrow = 1)
    mu_1 <- matrix(c(1.5, 2), nrow = 1)
    Sigma <- diag(1, nrow = 2)
}
num_digits <- 4
truths_1 <- round(truths[, 1], digits = num_digits)
truths_2 <- round(truths[, 2], digits = num_digits)
truths <- data.frame(j_1 = truths_1, j_2 = truths_2, type = truths$type)

## get the truths
if (length(args$vimp_measure) == 1) {
  truth <- subset(truths, type == args$vimp_measure)
} else {
  truth <- truths
}

## set up static args
V <- 5

## set up dynamic args (need to replicate these nreps_total/nreps_per_job times)
# ns <- c(100, 500, seq(1000, 8000, by = 1000)) # for large n simulations
ns <- c(100, 300, seq(500, 4000, 500)) # for small n simulations
js <- c(1, 2)
nreps_per_combo <- args$nreps_total/args$nreps_per_job
param_grid <- expand.grid(mc_id = 1:nreps_per_combo, j = js, n = ns)

## get dynamic args
current_dynamic_args <- param_grid[job_id, ]
if (length(args$type) > 1) {
  true_vals <- truth[, currend_dynamic_args$j]
} else {
  true_vals <- truth[current_dynamic_args$j]
}

## set up SuperLearner library
xgb_tune_params <- list(max_depth = 1, shrinkage = 0.01)
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params, detailed_names = TRUE, name_prefix = "xgb")
# gam_tune_params <- list(deg.gam = 1)
# gam_learners <- create.Learner("SL.gam", tune = gam_tune_params, detailed_names = TRUE, name_prefix = "gam")
learner_lib <- c(xgb_learners$names, "SL.glm", "SL.randomForest", "SL.mean")
## ---------------------------------------------
## run the simulation nreps_per_job times
## ---------------------------------------------
current_seed <- current_dynamic_args$mc_id + current_dynamic_args$j + ifelse(args$risk_type == "expected_loss", 0, 1) + current_dynamic_args$n + job_id
print(current_seed)
set.seed(current_seed)
system.time(sim_output <- replicate(args$nreps_per_job,
                        run_sim_binary_bivariate_once(n = current_dynamic_args$n,
                                                      j = current_dynamic_args$j,
                                                      p = p,
                                                      mu_0 = mu_0, mu_1 = mu_1, sigma = Sigma,
                                                      truth = true_vals,
                                                      b = args$b, V = V,
                                                      risk_type = args$risk_type,
                                                      learner_lib = learner_lib,
                                                      type = args$vimp_measure,
                                                      cv = args$cv),
                        simplify = FALSE))
## make a nice tibble, with mc_id
sim_output <- lapply(as.list(1:length(sim_output)), function(x) tibble::add_column(sim_output[[x]], mc_id = x))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0(args$sim_name, "_", paste0(args$vimp_measure, collapse = "_"), "_cv_", args$cv, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)
