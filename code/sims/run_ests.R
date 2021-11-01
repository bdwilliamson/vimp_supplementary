# run estimators with and without cross-fitting; no sample-splitting

# runs the estimation procedure based on a set of folds:
# if all folds are equal, then runs them on everything;
# otherwise, does cross-fitting
# @param draw the dataset
# @param learner the learner to fit
# @param reduced is this the reduced regression (based on a subset of the covariates)?
# @param folds the folds, for cross-fitting
# @return a list with the fitted values (on the validation fold) and a list of
#    predicted values
run_estimator_cv <- function(draw, learner = "gam", reduced = FALSE,
                             folds = rep(1, length(draw$y_cat)),
                             rf_opts = list(num_trees = 500, min_node_size = 3, 
                                            max_depth = NULL),
                             sl_opts = list(
                               fam = binomial(),
                               learners = c("SL.gam", "SL.ranger"),
                               method = "method.nnloglik"),
                             sample_splitting = FALSE, sample_splitting_folds = NULL) {
  this_x <- switch((reduced) + 1, draw$x, draw$red_x)
  if (sample_splitting) {
    V <- length(unique(folds)) / 2
    ss_V <- 2 * V
  } else {
    V <- length(unique(folds))
    ss_V <- V
  }
  if (learner == "SL") {
      # use CV.SuperLearner
      folds_lst <- lapply(as.list(sort(unique(folds))), function(v) {
        which(folds == v)
      })
      cv_sl_fit <- SuperLearner::CV.SuperLearner(
          Y = draw$y_cat, X = as.data.frame(this_x), SL.library = sl_opts$learners,
          family = sl_opts$fam, method = sl_opts$method,
          cvControl = list(V = ss_V, validRows = folds_lst),
          innerCvControl = list(list(V = 5))
      )
      fake_ss_fold <- switch((reduced) + 1, 1, 2)
      fit_lst <- vimp::extract_sampled_split_predictions(
        cvsl_obj = cv_sl_fit, sample_splitting = FALSE, full = !reduced,
        sample_splitting_folds = rep(fake_ss_fold, length(unique(folds)))
      )
      # if no sample splitting, extract all predictions
      if (!sample_splitting) {
          fit <- cv_sl_fit$SL.predict
          preds_lst <- vimp::extract_sampled_split_predictions(
              cvsl_obj = cv_sl_fit, sample_splitting = FALSE, full = !reduced,
              sample_splitting_folds = rep(fake_ss_fold, length(unique(folds)))
          )
      } else {
          fit <- cv_sl_fit$SL.predict
          if (is.null(sample_splitting_folds)) {
            sample_splitting_folds <- vimp::make_folds(
              y = unique(folds_lst), V = 2, C = rep(1, length(unique(folds_for_cv_sl)))
            )
          }
          preds_lst <- vimp::extract_sampled_split_predictions(
              cvsl_obj = cv_sl_fit, sample_splitting = TRUE,
              sample_splitting_folds = sample_splitting_folds, full = !reduced
          )
      }
  } else {
      # do cross-fitting manually
      fit <- vector("numeric", length = nrow(draw$x))
      fit_lst <- vector("list", length = ss_V)
      preds_lst <- vector("list", length = V)
      this_indx <- 1
      ss_fold <- switch((reduced) + 1, 1, 2)
      folds_lst <- vector("list", length = V)
      if (sample_splitting & is.null(sample_splitting_folds)) {
        sample_splitting_folds <- vimp::make_folds(y = unique(folds), V = 2,
                                                   C = rep(1, length(folds)))
      } else if (!is.null(sample_splitting_folds)) {
        # do nothing
      } else {
        sample_splitting_folds <- rep(1, ss_V)
      }
      for (v in seq_len(ss_V)) {
        train <- list(y_cat = draw$y_cat[folds != v],
                      x = draw$x[folds != v, ],
                      red_x = switch((!is.null(ncol(draw$red_x))) + 1,
                                     draw$red_x[folds != v],
                                     draw$red_x[folds != v, ]))
        test <- data.frame(y_cat = draw$y_cat[folds == v],
                           switch((!is.null(ncol(this_x))) + 1,
                                  this_x[folds == v],
                                  this_x[folds == v, ]))
        names(test) <- c("y", paste0("X", 1:(ncol(test) - 1)))
        fit_v <- run_estimator(train, test_draw = test,
                               learner = learner, reduced = reduced,
                               sl_opts = sl_opts, rf_opts = rf_opts)
        fit_lst[[v]] <- fit_v$fitted
        if (sample_splitting_folds[v] == ss_fold | !sample_splitting) {
          fit[folds == v] <- fit_v$fitted
          preds_lst[[this_indx]] <- fit_v$fitted
          folds_lst[[this_indx]] <- which(folds == v)
          this_indx <- this_indx + 1
        }
      }
  }
  list(fitted = fit, fitted_lst = fit_lst, preds_lst = preds_lst, 
       ss_folds = sample_splitting_folds, folds_lst = folds_lst)
}

# run the estimation procedure within a given fold
# @param draw the dataset
# @param test_draw the test dataset
# @param learner the learner to fit
# @param reduced is this the reduced regression (based on a subset of the covariates)?
# @return a list with the fitted object and a list of
#    predicted values
run_estimator <- function(draw, test_draw = NULL, learner = "gam",
                          reduced = FALSE, sl_opts = list(fam = binomial(), 
                                                          learners = c("SL.gam", 
                                                                       "SL.ranger", 
                                                                       method = "method.nnloglik")),
                          rf_opts = list(num_trees = 500, min_node_size = 3, 
                                         max_depth = NULL)) {
  this_x <- switch((reduced) + 1, draw$x, draw$red_x)
  dat <- data.frame(y = draw$y_cat, this_x)
  names(dat) <- c("y", paste0("X", 1:(ncol(dat) - 1)))
  if (learner == "glm") {
    fit <- glm(y ~ ., family = sl_opts$fam, data = dat)
    if (is.null(test_draw)) {
      preds <- as.numeric(predict(fit, type = "response"))
    } else {
      preds <- as.numeric(predict(fit, newdata = test_draw, type = "response"))
    }
    coefs <- coef(fit)
    fit_out <- paste0(names(coefs), ": ", round(coefs, 3), collapse = "; ")
  } else if (learner == "gam") {
    x_nms <- names(dat)[-1]
    gam_form <- as.formula(paste0("y ~ ", paste0("s(", x_nms, ", bs = 'tp')",
                                                 collapse = "+")))
    fit <- mgcv::gam(gam_form, family = sl_opts$fam, data = dat,
                     method = "REML",
                     optimizer = c("outer", "bfgs"))
    if (is.null(test_draw)) {
      preds <- as.numeric(predict(fit, type = "response"))
    } else {
      preds <- as.numeric(predict(fit, newdata = test_draw, type = "response"))
    }
    coefs <- coef(fit)
    fit_out <- paste0(names(coefs), ": ", round(coefs, 3), collapse = "; ")
  } else if (learner == "ranger") {
    newy <- as.factor(dat$y)
    newdat <- data.frame(Y = newy, dat[, -1, drop = FALSE])
    fit <- ranger::ranger(Y ~ ., data = newdat, num.trees = rf_opts$num_trees,
                          replace = TRUE, probability = (sl_opts$fam$family == "binomial"),
                          num.threads = 1, sample.fraction = 1,
                          min.node.size = rf_opts$min_node_size, 
                          max.depth = rf_opts$max_depth,
                          mtry = switch(
                            ((ncol(dat) - 1) == 1) + 1, floor(sqrt(ncol(dat))),
                            ncol(dat) - 1
                          ))
    if (is.null(test_draw)) {
      preds <- predict(fit, data = dat[, -1, drop = FALSE])$predictions[, "1"]
    } else {
      preds <- predict(fit, data = test_draw[, -1, drop = FALSE])$predictions[, "1"]
    }
    fit_out <- paste0("ntrees: ", rf_opts$num_trees, "; min.node.size: ", 
                      rf_opts$min_node_size, "; max.depth: ", rf_opts$max_depth)
  } else if (learner == "xgboost") {
    x_mat <- model.matrix(~ . - 1, dat[, -1, drop = FALSE])
    xgb_mat <- xgboost::xgb.DMatrix(data = x_mat, label = dat$y)
    fit <- xgboost::xgboost(data = xgb_mat, 
                            objective = switch((sl_opts$fam$family == "binomial") + 1, 
                                               "reg:squarederror",
                                               "binary:logistic"),
                            params = list(eta = rf_opts$eta, 
                                          max_depth = rf_opts$max_depth),
                            nrounds = 1000, verbose = 0)
    if (is.null(test_draw)) {
      preds <- predict(fit, newdata = x_mat)
    } else {
      new_mat <- model.matrix(~ . - 1, test_draw[, -1, drop = FALSE])
      preds <- predict(fit, newdata = new_mat)
    }
    fit_out <- paste0("shrinkage: ", rf_opts$eta, "; max_depth: ", rf_opts$max_depth)
  } else if (learner == "SL") {
    x <- dat[, -1, drop = FALSE]
    fit <- SuperLearner::SuperLearner(Y = dat$y, X = x, cvControl = list(V = 5),
                                      family = sl_opts$fam, method = sl_opts$method,
                                      SL.library = sl_opts$learners)
    if (is.null(test_draw)) {
      preds <- fit$SL.predict
    } else {
      preds <- SuperLearner::predict.SuperLearner(
        object = fit, newdata = test_draw[, -1, drop = FALSE], onlySL = TRUE
      )$pred
    }
    out_df <- rbind.data.frame(round(fit$cvRisk, 3), round(fit$coef, 3))
    names(out_df) <- names(fit$coef)
    fit_out <- paste0(names(out_df), ": ", apply(out_df, 2, function(x) paste0(x, collapse = "/")),
                      collapse = "; ")
} else {
    stop(paste0("The requested learner is not currently implemented. ",
    "Please enter one of 'glm', 'gam', 'ranger', or 'SL'."))
  }
  list(fit = fit, fitted = preds, fit_out = fit_out)
}
