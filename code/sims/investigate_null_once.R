# run through a simple estimation procedure one time
do_one <- function(iteration = 1, n = 100, p = 2, j = 1, x_mean = rep(0, p),
                   beta_0 = matrix(rep(0, p)), Sigma = diag(1, p), truth = 0,
                   cv = FALSE, learner = "gam", sl_lib = c("SL.gam", "SL.ranger"),
                   type = "auc", V = 5, bootstrap = FALSE, linkfun = "logit",
                   b = 1000) {
  # generate data, make folds
  draw <- gen_probit_data(nsim = n, p = p, x_mean = x_mean, Sigma = Sigma,
                          beta_0 = beta_0, j = j)
  # note that for cross-fitting, we use 2V folds (though only evaluate on
  # V folds for both full and reduced)
  cv_folds <- vimp::make_folds(y = draw$y_cat, V = 2 * V,
                               C = rep(1, length(draw$y_cat)), stratified = TRUE)
  ss_folds <- vimp::make_folds(y = unique(cv_folds), V = 2,
                               C = rep(1, length(unique(cv_folds))))
  sl_opts <- get_sl_opts(y = draw$y_cat, link_fun = linkfun, learners = sl_lib)
  # if learner == "rf", use CV to get optimal params
  if (learner == "ranger") {
    num_trees <- c(500, 1000, 1500, 2000, 5000)
    min_node_size <- c(10)
    max_depth <- c(1, 3, 5)
    rf_opts_full <- get_rf_opts(draw = draw, family = sl_opts$fam,
                           num_trees = num_trees, min_node_size = min_node_size,
                           max_depth = max_depth)
    rf_opts_redu <- get_rf_opts(draw = list(y_cat = draw$y_cat,
                                            x = draw$red_x),
                                family = sl_opts$fam,
                                num_trees = num_trees, min_node_size = min_node_size,
                                max_depth = max_depth)
  } else if (learner == "xgboost") {
    rf_opts_full <- get_xgb_opts(draw = draw, family = sl_opts$fam, shrinkage = 0.1,
                            max_depth = c(2, 4))
    rf_opts_redu <- get_xgb_opts(draw = draw, family = sl_opts$fam, shrinkage = 0.1,
                                 max_depth = c(2, 4))
  } else {
    rf_opts_full <- rf_opts_redu <- NULL
  }
  # obtain estimators of the full and reduced regression functions
  # if cv, do this using cross-fitting
  if (cv) {
    cv_full_fit_lst <- run_estimator_cv(draw, learner = learner, reduced = FALSE,
                                     folds = cv_folds, sl_opts = sl_opts,
                                     rf_opts = rf_opts_full,
                                     sample_splitting = TRUE,
                                     sample_splitting_folds = ss_folds)
    cv_redu_fit_lst <- run_estimator_cv(draw, learner = learner, reduced = TRUE,
                                     folds = cv_folds, sl_opts = sl_opts,
                                     rf_opts = rf_opts_redu,
                                     sample_splitting = TRUE,
                                     sample_splitting_folds = ss_folds)
    cv_full_fits_lst <- cv_full_fit_lst$fitted_lst
    cv_redu_fits_lst <- cv_redu_fit_lst$fitted_lst
    cv_full_preds_lst <- cv_full_fit_lst$preds_lst
    cv_redu_preds_lst <- cv_redu_fit_lst$preds_lst
  }
  # fitted values for non-cross-fitted SE
  cf_folds_1 <- sort(unique(cv_folds))[ss_folds == 1]
  full_data_ss_folds <- ifelse(cv_folds %in% cf_folds_1, 1, 2)
  draw_1 <- get_sample_split_data(draw = draw, ss_folds = full_data_ss_folds, 
                                  full = TRUE)
  full_fit_lst <- run_estimator(draw = draw_1, test_draw = NULL,
                                learner = learner, reduced = FALSE,
                                rf_opts = rf_opts_full, sl_opts = sl_opts)
  draw_2 <- get_sample_split_data(draw = draw, ss_folds = full_data_ss_folds, 
                                  full = FALSE)
  redu_fit_lst <- run_estimator(draw = draw_2, test_draw = NULL,
                                learner = learner, reduced = TRUE,
                                rf_opts = rf_opts_redu, sl_opts = sl_opts)
  full_preds <- full_fit_lst$fitted
  redu_preds <- redu_fit_lst$fitted

  if (cv) {
    ests_cf_se <- lapply(as.list(type), function(t) {
      cv_vim(Y = draw$y_cat, X = draw$x, cross_fitted_f1 = cv_full_preds_lst,
             cross_fitted_f2 = cv_redu_preds_lst, indx = j, V = V, 
             sample_splitting = TRUE, cross_fitting_folds = cv_folds,
             type = t, run_regression = FALSE, alpha = 0.05,
             scale = "identity", delta = 0, sample_splitting_folds = ss_folds,
             scale_est = FALSE)
    })
    ests_non_cf_se <- lapply(as.list(type), function(t) {
      cv_vim(Y = draw$y_cat, X = draw$x, cross_fitted_f1 = cv_full_preds_lst,
             cross_fitted_f2 = cv_redu_preds_lst, f1 = full_preds, f2 = redu_preds,
             indx = j, V = V, sample_splitting = TRUE, cross_fitting_folds = cv_folds,
             type = t, run_regression = FALSE, alpha = 0.05,
             scale = "identity", delta = 0, sample_splitting_folds = ss_folds,
             scale_est = FALSE, cross_fitted_se = FALSE)
    })
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
          scale = "identity", delta = 0, sample_splitting = TRUE,
          sample_splitting_folds = full_data_ss_folds, scale_est = FALSE)
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
