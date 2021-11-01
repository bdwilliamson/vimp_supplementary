#!/usr/local/bin/Rscript

# a simple simulation to see study the properties of our proposal under the
# alternative (both features are important)

#### load required libraries and functions ####
library("methods")
lib_paths <- .libPaths()
library("argparse")
library("ranger")
library("xgboost")
library("dplyr")
library("tibble")
library("data.table")
library("SuperLearner")
library("nloptr")
library("here")

# edit the following line if using a different system than Slurm
job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(job_id)) { # if running locally
  job_id <- 1
  # require vimp version >= 2.2.3
  library("vimp")
  code_prefix <- "sim_code"
  prefix <- "sim_output/"
} else {
  # require vimp version >= 2.2.3
  library("vimp", lib.loc = switch((length(lib_paths) > 1) + 1, NULL, lib_paths[2]))
  code_prefix <- "."
  # edit the following line to where you wish to save results
  prefix <- paste0("path_to_results")
}

source(here(code_prefix, "gen_data.R"))
source(here(code_prefix, "investigate_null_once.R"))
source(here(code_prefix, "run_ests.R"))
source(here(code_prefix, "utils.R"))

#### pull in command-line arguments, set up the simulation ####
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "investigate-null-props",
                    help = "the name of the simulation")
parser$add_argument("--learner", default = "glm",
                    help = "which estimator to run (glm, gam, ranger, SL)")
parser$add_argument("--nreps-total", type = "double", default = 500,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 1,
                    help = "number of replicates for each job")
parser$add_argument("--cv", type = "double", default = 0,
                    help = "run cross-validated estimators?")
parser$add_argument("--j", default = 3, type = "double",
                    help = "which variable to estimate importance for")
parser$add_argument("--p", default = 4, type = "double",
                    help = "dimension to consider")
parser$add_argument("--corr", default = 0, type = "double",
                    help = "should we include correlation between predictors?")
parser$add_argument("--bootstrap", default = 0, type = "double",
                    help = "should we use the bootstrap for variance estimation?")
args <- parser$parse_args()
print(args)

if (grepl("SL", args$learner)) {
  # use gam gam
  library("gam")
  if (args$p > 2) {
    library("glmnet")
  }
} else {
  # use mgcv gam
  library("mgcv")
}

#### true parameter values ####
x_mean <- rep(0, args$p)
beta_0 <- matrix(c(2.5, 3.5, rep(0, args$p - 2)))
Sigma <- diag(1, args$p)
# covariances; only used if args$corr > 0
rho_13 <- 0.7
rho_24 <- 0.2
if (as.logical(args$corr)) {
  Sigma[1, 3] <- Sigma[3, 1] <- rho_13
  Sigma[2, 4] <- Sigma[4, 2] <- rho_24
}

truth <- get_true_values(args, num_digits = 4, dir = code_prefix)

#### set up static, dynamic args ####
V <- 5
# sample sizes and features considered depend on the dimension
if (args$p > 4) {
  ns <- c(500, 3000)
  if (args$j == 24) {
    args$j <- "2,4"
  } else if (args$j == 13) {
    args$j <- "1,3"
  }
} else {
  ns <- seq(250, 2000, 250)
}
nreps_per_combo <- args$nreps_total / args$nreps_per_job
param_grid <- expand.grid(mc_id = 1:nreps_per_combo, j = args$j, n = ns)

# get dynamic args
current_dynamic_args <- param_grid[job_id, ]
if (current_dynamic_args$j == "2,4") {
  true_vals <- truth$j_24_group
} else if (current_dynamic_args$j == "1,3") {
  true_vals <- truth$j_13_group
} else {
  true_vals <- truth %>%
    dplyr::pull(paste0("j_", current_dynamic_args$j))
}

# set up a SL library (only used if learner == "SL")
if (args$p <= 4) {
  max_depth <- 4
} else {
  max_depth <- 1
}
xgb_tune_params <- list(max_depth = max_depth, shrinkage = 0.1)
xgb_learners <- create.Learner("SL.xgboost", tune = xgb_tune_params,
                               detailed_names = TRUE, name_prefix = "xgb")
# only use glmnet if p > 2
learner_lib <- c(xgb_learners$names, "SL.gam", "SL.ranger", "SL.glm")
if (args$p > 2 & args$p <= 4) {
  learner_lib <- c(learner_lib, "SL.glmnet")
} else {
  learner_lib <- c(xgb_learners$names, "SL.ranger", "SL.glmnet")
}

# allow for group importance
if (is.factor(current_dynamic_args$j)) {
  this_j <- as.numeric(unlist(strsplit(as.character(current_dynamic_args$j),
                                       ",", fixed = TRUE)))
} else {
  this_j <- current_dynamic_args$j
}

#### run the simulation nreps_per_job times, save output ####
if (args$p > 4) {
  j_seed_contrib <- case_when(
    args$j == "1" ~ 0,
    args$j == "2" ~ 1,
    args$j == "3" ~ 2,
    args$j == "4" ~ 3,
    args$j == "1,3" ~ 4,
    args$j == "2,4" ~ 5,
    args$j == "5" ~ 6
  ) * 1000
} else {
  j_seed_contrib <- switch((args$j == 2) + 1,
                           as.numeric(current_dynamic_args$j) * 1000,
                           as.numeric(current_dynamic_args$j))
}
current_seed <- current_dynamic_args$mc_id +
  j_seed_contrib +
  current_dynamic_args$n * 1000

print(current_seed)
set.seed(current_seed)
system.time(sim_output <- sapply(
  1:args$nreps_per_job,
  function(i) {
    do_one(iteration = i + args$nreps_per_job * (current_dynamic_args$mc_id - 1),
           n = current_dynamic_args$n, p = args$p, j = this_j,
           x_mean = x_mean, beta_0 = beta_0, Sigma = Sigma,
           type = c("accuracy", "auc"), truth = true_vals,
           cv = as.logical(args$cv), learner = args$learner,
           linkfun = "probit",
           sl_lib = learner_lib, V = V, bootstrap = as.logical(args$bootstrap))
  }, simplify = FALSE)
)
# make a nice tibble, with mc_id
sim_output_tib <- tibble::as_tibble(rbindlist(sim_output))
# note that if SE = 0, test should be FALSE
sim_output_tib %>%
  mutate(bias_init = sqrt(n) * (est - truth),
         cover_init = cil <= truth & ciu >= truth,
         test_fixed = ifelse(is.na(test & se == 0), FALSE, test)) %>%
  group_by(type, j) %>%
  summarize(mn_est = mean(est), bias = mean(bias_init),
            cover = mean(cover_init),
            type_1_error = mean(test_fixed))
file_name <- paste0(args$learner, "_cv_", args$cv, "_j_",
                    args$j, "_p_", args$p, "_corr_", args$corr,
                    "_boot_", args$bootstrap, "_",
                    job_id, ".rds")
output_dir <- paste0(prefix, args$sim_name, "/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
saveRDS(sim_output_tib, file = paste0(output_dir, file_name))
print("Simulation finished!")
