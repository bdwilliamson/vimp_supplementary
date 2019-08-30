#!/usr/local/bin/Rscript
## ------------------------------------------------
## FILE: run_sim_binary_bivariate_once.R
## CREATED: 06 March 2019 by Brian Williamson
## PURPOSE: run the simulation once:
##          (a) create a dataset,
##          (b) run the estimators,
##          (c) return results
## ------------------------------------------------
run_sim_binary_bivariate_once <- function(n, j, p,
                                          mu_0, mu_1, sigma, 
                                          truth,
                                          b, V, risk_type, learner_lib,
                                          type, cv) {
    ## generate data
    draw <- gen_data(a = n, mu_0 = mu_0, mu_1 = mu_1, sigma = sigma, p = p, j = j)
    ## run the estimators on the data
    est <- nonparametric_estimator(draw = draw, V = V, learner_lib = learner_lib, 
                                   type = type, cv = cv)
    ## make a tibble with the output; names are based on type
    out <- tibble::tibble(n = n, j = j, truth = unlist(truth))
    out <- tibble::add_column(out, est = unlist(lapply(est, function(x) x$est)))
    out <- tibble::add_column(out, se = unlist(lapply(est, function(x) x$se)))                      
    out <- tibble::add_column(out, cil = unlist(lapply(est, function(x) x$ci[1])))
    out <- tibble::add_column(out, ciu = unlist(lapply(est, function(x) x$ci[2])))
    out <- tibble::add_column(out, test = unlist(lapply(est, function(x) x$test)))
    out <- tibble::add_column(out, p_value = unlist(lapply(est, function(x) x$p_value)))
    out <- tibble::add_column(out, risk_full = unlist(lapply(est, function(x) x$predictiveness_full)))
    out <- tibble::add_column(out, risk_reduced = unlist(lapply(est, function(x) x$predictiveness_reduced)))
    out <- tibble::add_column(out, type = type)
    
    ## if risk_type != "expected_loss", then also run the naive estimator
    if (risk_type != "expected_loss") {
      naive <- nonparametric_naive_deviance(draw, V, learner_lib, cv = FALSE)
      boot_df <- tibble(y = draw$y, x = draw$x, y_cat = draw$y_cat, red_x = draw$red_x)
      naives <- boot::boot(boot_df, nonparametric_naive_deviance_boot, R = b, V = V, learner_lib = learner_lib, cv = FALSE)
      naive_ci <- boot::boot.ci(naives, conf = 0.95, type = "perc")$percent[4:5]
      out <- tibble::add_column(out, naive_deviance = naive, naive_cil_deviance = naive_ci[1], naive_ciu_deviance = naive_ci[2])
    }
    
    return(out)
}