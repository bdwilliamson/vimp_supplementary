#!/usr/local/bin/Rscript
## ------------------------------------------------
## FILE: sim_no_ss_ests.R
## CREATED: 31 March 2020 by Brian Williamson
## PURPOSE: VIM estimates for no ss sim
## ------------------------------------------------

## ----------------------------------------------------------------------------------------
## nonparametric estimators
## ----------------------------------------------------------------------------------------
estimator <- function(draw, V, learner_lib, type, cv, delta) {
  if (class(draw$x) != "data.frame") {
    draw$x <- as.data.frame(draw$x)
    draw$red_x <- as.data.frame(draw$red_x)
    names(draw$red_x) <- "V1"
  } else {
    draw$red_x <- as.data.frame(draw$red_x)
    names(draw$red_x) <- "V1"
  }
  if (cv) {
    ## set up the cv folds
    inner_folds_1 <- make_folds(draw, V = V, stratified = TRUE)
    inner_folds_2 <- make_folds(draw, V = V, stratified = TRUE)

    ## initialize lists for the fitted values, run the regressions
    full_fit <- list()
    redu_fit <- list()
    for (v in 1:V) {
      full_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[inner_folds_1 != v], X = draw$x[inner_folds_1 != v, , drop = FALSE], family = binomial(), cvControl = list(V = 10), SL.library = learner_lib, method = "method.CC_nloglik"))
      if (any(is.na(full_mod$cvRisk))) {
        full_fit[[v]] <- NA
      } else {
        full_fit[[v]] <- SuperLearner::predict.SuperLearner(full_mod, newdata = draw$x[inner_folds_1 == v, , drop = FALSE], onlySL = TRUE)$pred[, 1]
      }
      redu_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[inner_folds_2 != v], X = draw$red_x[inner_folds_2 != v, , drop = FALSE], family = binomial(), cvControl = list(V = 10), SL.library = learner_lib, method = "method.CC_nloglik")) # remove "method.NNloglik" due to numerical errors
      if (any(is.na(redu_mod$cvRisk))) {
        redu_fit[[v]] <- NA
      } else {
        redu_fit[[v]] <- SuperLearner::predict.SuperLearner(redu_mod, newdata = draw$red_x[inner_folds_2 == v, , drop = FALSE], onlySL = TRUE)$pred[, 1]
      }
    }
    ests <- lapply(as.list(type), function(x) est <- get_est_ci_test(y = draw$y_cat, x = draw$x, f1 = full_fit, f2 = redu_fit, alpha = 0.05, delta = delta, type = x, cv = TRUE, folds = list(folds_1 = inner_folds_1, folds_2 = inner_folds_2)))
  } else {
    draw_1 <- draw
    draw_2 <- draw
    ## estimate the regression functions
    full_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw_1$y_cat, X = draw_1$x, family = binomial(), cvControl = list(V = 10), SL.library = learner_lib, method = "method.CC_nloglik"))
    full_fit <- full_mod$SL.predict
    full <- cbind(1 - full_fit, full_fit) # predicting Y = 1
    redu_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw_2$y_cat, X = draw_2$red_x, family = binomial(), cvControl = list(V = 10), SL.library = learner_lib, method = "method.CC_nloglik")) # remove "method.NNloglik" due to numerical errors
    redu_fit <- redu_mod$SL.predict
    redu <- cbind(1 - redu_fit, redu_fit)

    ests <- lapply(as.list(type), function(x) est <- get_est_ci_test(y = draw$y_cat, x = draw$x, f1 = full_fit, f2 = redu_fit, alpha = 0.05, delta = delta, type = x, cv = FALSE))
  }

  return(ests)
}

make_folds <- function(draw, V, stratified = TRUE) {
  if (stratified) {
    y_1 <- draw$y_cat == 1
    y_0 <- draw$y_cat == 0
    folds_1 <- rep(seq_len(V), length = sum(y_1))
    folds_1 <- sample(folds_1)
    folds_0 <- rep(seq_len(V), length = sum(y_0))
    folds_0 <- sample(folds_0)
    folds <- vector("numeric", length(draw$y_cat))
    folds[y_1] <- folds_1
    folds[y_0] <- folds_0
  } else {
    folds <- rep(seq_len(V), length = length(draw$y_cat))
    folds <- sample(folds)
  }
  return(folds)
}

## obtain estimates, CIs, hypothesis test of 0.05 null for a given type of VIM
get_est_ci_test <- function(y, x, f1, f2, alpha, delta, type, cv, folds) {
  if (cv) {
    point_est_full <- cv_predictiveness_point_est(fitted_values = f1, y = y, folds = folds[[1]], type = type)$point_est
    point_est_redu <- cv_predictiveness_point_est(fitted_values = f2, y = y, folds = folds[[2]], type = type)$point_est
    ic_full <- cv_predictiveness_update(fitted_values = f1, y = y, folds = folds[[1]], type = type)$ic
    ic_redu <- cv_predictiveness_update(fitted_values = f2, y = y, folds = folds[[2]], type = type)$ic
  } else {
    point_est_full <- predictiveness_point_est(fitted_values = f1, y = y, type = type)
    point_est_redu <- predictiveness_point_est(fitted_values = f2, y = y, type = type)
    ic_full <- predictiveness_update(fitted_values = f1, y = y, type = type)
    ic_redu <- predictiveness_update(fitted_values = f2, y = y, type = type)
  }
  se_full <- sqrt(mean(ic_full ^ 2))
  se_redu <- sqrt(mean(ic_redu ^ 2))
  vimp_est <- point_est_full - point_est_redu
  vimp_est_ic <- ic_full - ic_redu
  vimp_est_se <- vimp_se(vimp_est, vimp_est_ic, scale = "identity")
  vimp_est_ci <- vimp_ci(vimp_est, vimp_est_se, scale = "identity", level = 1 - alpha)
  vimp_tstat <- (vimp_est - delta) * (se_full ^ 2 / length(ic_full) + se_redu ^ 2 / length(ic_redu)) ^ -(1/2)
  vimp_pval <- 1 - pnorm(vimp_tstat)
  vimp_test <- vimp_pval < alpha
  return(list(est = vimp_est, se = vimp_est_se, ci = vimp_est_ci, p_value = vimp_pval, test = vimp_test,
              point_est_full = point_est_full, point_est_redu = point_est_redu))
}
