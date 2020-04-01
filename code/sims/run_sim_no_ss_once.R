#!/usr/local/bin/Rscript
## ------------------------------------------------
## FILE: run_sim_no_ss_once.R
## CREATED: 31 March 2020 by Brian Williamson
## PURPOSE: run the simulation once:
##          (a) create a dataset,
##          (b) run the estimators,
##          (c) return results
## ------------------------------------------------
run_sim_no_ss_once <- function(n, j, p, mu_0, mu_1, sigma, truth,
                                          b, V, risk_type, learner_lib,
                                          type, cv, delta = 0.05) {
    ## generate data
    draw <- gen_data(a = n, mu_0 = mu_0, mu_1 = mu_1, sigma = sigma, p = p, j = j)
    ## run the estimators on the data
    est_lst <- estimator(draw = draw, V = V, learner_lib = learner_lib, type = type, cv = cv, delta = delta)
    ## make a tibble with the output; names are based on type
    out <- tibble::tibble(n = n, j = j, truth = unlist(truth))
    out <- tibble::add_column(out, est = unlist(lapply(est_lst, function(x) x$est)),
                              se = unlist(lapply(est_lst, function(x) x$se)),
                              cil = unlist(lapply(est_lst, function(x) x$ci[1])),
                              ciu = unlist(lapply(est_lst, function(x) x$ci[2])),
                              test = unlist(lapply(est_lst, function(x) x$test)),
                              p_value = unlist(lapply(est_lst, function(x) x$p_value)),
                              risk_full = unlist(lapply(est_lst, function(x) x$point_est_full)),
                              risk_reduced = unlist(lapply(est_lst, function(x) x$point_est_redu)),
                              type = type, delta = delta)
    return(out)
}
