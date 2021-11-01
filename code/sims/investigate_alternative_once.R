# run through a simple estimation procedure one time
do_one <- function(iteration = 1, n = 100, p = 2, j = 1, x_mean = rep(0, p),
                   beta_0 = matrix(rep(0, p)), Sigma = diag(1, p), truth = 0,
                   cv = FALSE, learner = "gam", sl_lib = c("SL.gam", "SL.ranger"),
                   type = "auc", V = 5, bootstrap = FALSE, linkfun = "logit",
                   b = 1000, boot_type = "perc") {
  # generate data, make folds
  draw <- gen_probit_data(nsim = n, p = p, x_mean = x_mean, Sigma = Sigma,
                          beta_0 = beta_0, j = j)
  # note that for cross-fitting, we use V folds
  # since we're not doing sample-splitting
  cv_folds <- vimp::make_folds(draw$y_cat, V = V,
                               C = rep(1, length(draw$y_cat)), stratified = TRUE)
  sl_opts <- get_sl_opts(y = draw$y_cat, link_fun = linkfun, learners = sl_lib)
  # if learner == "rf", use CV to get optimal params
  if (learner == "ranger") {
    rf_opts <- get_rf_opts(draw = draw, family = sl_opts$fam,
                           num_trees = c(100, 500, 1000),
                           min_node_size = 1,
                           max_depth = 5)
  } else if (learner == "xgboost") {
    rf_opts <- get_xgb_opts(draw = draw, family = sl_opts$fam, shrinkage = 0.1,
                            max_depth = c(2, 4))
  } else {
    rf_opts <- NULL
  }
  # obtain estimators of the full and reduced regression functions
  # if cv, do this using cross-fitting
  if (cv) {
    cv_full_fit_lst <- run_estimator_cv(draw, learner = learner, reduced = FALSE,
                                     folds = cv_folds, sl_opts = sl_opts,
                                     rf_opts = rf_opts,
                                     sample_splitting = FALSE)
    cv_redu_fit_lst <- run_estimator_cv(draw, learner = learner, reduced = TRUE,
                                     folds = cv_folds, sl_opts = sl_opts,
                                     rf_opts = rf_opts,
                                     sample_splitting = FALSE)
    cv_full_fits_lst <- cv_full_fit_lst$fitted_lst
    cv_redu_fits_lst <- cv_redu_fit_lst$fitted_lst
    cv_full_preds_lst <- cv_full_fit_lst$preds_lst
    cv_redu_preds_lst <- cv_redu_fit_lst$preds_lst
  }
  # regardless, fit to the whole dataset
  full_fit_lst <- run_estimator(draw = draw, learner = learner, reduced = FALSE,
                                rf_opts = rf_opts, sl_opts = sl_opts)
  redu_fit_lst <- run_estimator(draw = draw, learner = learner, reduced = TRUE,
                                rf_opts = rf_opts, sl_opts = sl_opts)
  full_preds <- full_fit_lst$fitted
  redu_preds <- redu_fit_lst$fitted

  if (cv) {
    if (!bootstrap) {
      ests_cf_se <- lapply(as.list(type), function(t) {
        cv_vim(Y = draw$y_cat, X = draw$x, cross_fitted_f1 = cv_full_preds_lst,
               cross_fitted_f2 = cv_redu_preds_lst, f1 = cv_full_fits_lst, f2 = cv_redu_fits_lst,
               indx = j, V = V, sample_splitting = FALSE, cross_fitting_folds = cv_folds,
               type = t, run_regression = FALSE, alpha = 0.05,
               scale = "identity", delta = 0, bootstrap = bootstrap, b = b,
               cross_fitted_se = ifelse(bootstrap, FALSE, TRUE))
      })
      ests_non_cf_se <- lapply(as.list(type), function(t) {
        cv_vim(Y = draw$y_cat, X = draw$x, cross_fitted_f1 = cv_full_preds_lst,
               cross_fitted_f2 = cv_redu_preds_lst, f1 = full_preds, f2 = redu_preds,
               indx = j, V = V, sample_splitting = FALSE, cross_fitting_folds = cv_folds,
               type = t, run_regression = FALSE, alpha = 0.05,
               scale = "identity", delta = 0, bootstrap = bootstrap, b = b,
               cross_fitted_se = FALSE)
      })  
    } else {
      ests_cf_se <- lapply(as.list(type), function(t) {
        cv_vim(Y = draw$y_cat, X = draw$x, cross_fitted_f1 = cv_full_preds_lst,
               cross_fitted_f2 = cv_redu_preds_lst, f1 = full_preds, f2 = redu_preds,
               indx = j, V = V, sample_splitting = FALSE, cross_fitting_folds = cv_folds,
               type = t, run_regression = FALSE, alpha = 0.05,
               scale = "identity", delta = 0, bootstrap = bootstrap, b = b,
               cross_fitted_se = FALSE, boot_interval_type = boot_type)
      })
      ests_non_cf_se <- ests_cf_se
    }
    est_vim <- data.table::rbindlist(lapply(ests_cf_se, function(l) l$mat))
    est_full_pred <- do.call(rbind, lapply(ests_cf_se, function(l) l$predictiveness_full))
    est_redu_pred <- do.call(rbind, lapply(ests_cf_se, function(l) l$predictiveness_reduced))
    non_cf_se <- do.call(c, lapply(ests_non_cf_se, function(l) l$se))
    est_vim <- est_vim %>%
      dplyr::bind_cols(non_cf_se = non_cf_se)
  } else {
    ests <- lapply(as.list(type), function(t) {
      vim(Y = draw$y_cat, X = draw$x, f1 = full_preds, f2 = redu_preds,
          indx = j, type = t, run_regression = FALSE, alpha = 0.05,
          scale = "identity", delta = 0, sample_splitting = FALSE,
          bootstrap = bootstrap, b = b)
    })
    est_vim <- data.table::rbindlist(lapply(ests, function(l) l$mat))
    est_full_pred <- do.call(rbind, lapply(ests, function(l) l$predictiveness_full))
    est_redu_pred <- do.call(rbind, lapply(ests, function(l) l$predictiveness_reduced))
    est_vim$non_cf_se <- est_vim$se
  }
  # return output
  if (p > 2) {
    corr <- 1 - as.numeric(Sigma[1, 3] == 0)
  } else {
    corr <- 0
  }
  dplyr::bind_cols(tibble::tibble(mc_id = iteration, n = n, p = p, corr = corr,
                                  j = paste0(j, collapse = ","), estimator = learner,
                                  cv = as.numeric(cv), type = type, truth = truth),
                   est_vim %>% select(-s),
                   tibble::tibble(full_pred = as.numeric(est_full_pred),
                                  reduced_pred = as.numeric(est_redu_pred),
                                  fit_out = full_fit_lst$fit_out))
}
