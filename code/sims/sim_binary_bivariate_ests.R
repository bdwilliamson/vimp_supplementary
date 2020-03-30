#!/usr/local/bin/Rscript
## ------------------------------------------------
## FILE: sim_binary_bivariate_ests.R
## CREATED: 06 March 2019 by Brian Williamson
## PURPOSE: estimators of regression functions
##          for binary outcome, bivariate predictor
## ------------------------------------------------

## ----------------------------------------------------------------------------------------
## nonparametric estimators
## ----------------------------------------------------------------------------------------
nonparametric_estimator <- function(draw, V, learner_lib, type, cv) {
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
    folds <- rep(seq_len(V), length = length(draw$y_cat))
    folds <- sample(folds)
    
    ## initialize lists for the fitted values, run the regressions
    full_fit <- list()
    redu_fit <- list()
    for (v in 1:V) {
      full_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[folds != v], X = draw$x[folds != v, , drop = FALSE], family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik"))
      if (any(is.na(full_mod$cvRisk))) {
        full_fit[[v]] <- NA
      } else {
        full_fit[[v]] <- SuperLearner::predict.SuperLearner(full_mod, newdata = draw$x[folds == v, , drop = FALSE], onlySL = TRUE)$pred[, 1]  
      }
      
      
      redu_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[folds != v], X = draw$red_x[folds != v, , drop = FALSE], family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik")) # remove "method.NNloglik" due to numerical errors
      if (any(is.na(redu_mod$cvRisk))) {
        redu_fit[[v]] <- NA  
      } else {
        redu_fit[[v]] <- SuperLearner::predict.SuperLearner(redu_mod, newdata = draw$red_x[folds == v, , drop = FALSE], onlySL = TRUE)$pred[, 1]
      }
      
    }
    
    ests <- lapply(as.list(type), function(x) cv_vim(Y = draw$y_cat, X = draw$x, 
                                                       f1 = full_fit, f2 = redu_fit, 
                                                       indx = draw$j, V = V, folds = folds, 
                                                       type = x, run_regression = FALSE,
                                                       alpha = 0.05, na.rm = TRUE, scale = "identity"))
  } else {
    full_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat, X = draw$x, family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik"))
    full_fit <- full_mod$SL.predict
    full <- cbind(1 - full_fit, full_fit) # predicting Y = 1
    red_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat, X = draw$red_x, family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik")) # remove "method.NNloglik" due to numerical errors
    red_fit <- red_mod$SL.predict
    red <- cbind(1 - red_fit, red_fit)
    
    ## create splits and run on splits, for hypothesis testing
    folds <- rep(seq_len(2), length = length(full_fit))
    folds <- sample(folds)
    full_mod_split <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[folds == 1], X = draw$x[folds == 1, ], family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik"))
    redu_mod_split <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[folds == 2], X = draw$red_x[folds == 2, , drop = FALSE], family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik"))
    
    full_fit_split <- full_mod_split$SL.predict
    redu_fit_split <- redu_mod_split$SL.predict
    
    ests <- lapply(as.list(type), function(x) est <- vim(Y = draw$y_cat, X = draw$x, 
                                                           f1 = full_fit, f2 = red_fit, 
                                                           indx = draw$j, type = x, run_regression = FALSE,
                                                           alpha = 0.05, na.rm = TRUE, f1_split = full_fit_split,
                                                           f2_split = redu_fit_split, folds = folds, scale = "identity") )
  }

  return(ests)
}

## this function is only for the deviance (hard to imagine how to do this for the other parameters)
## here, rather than using an expected loss function as our risk,
## we use something slightly different (which blows things up); this highlights the 
## need to be careful about how you define your risk
nonparametric_naive_deviance <- function(draw, V, learner_lib, cv = TRUE) {
  ## fit (discrete, if type == "discrete_sl") super learner
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
    folds <- rep(seq_len(V), length = length(draw$y_cat))
    folds <- sample(folds)
    
    ## initialize lists for the fitted values, run the regressions
    full_fit <- list()
    redu_fit <- list()
    numerators <- vector("numeric", V)
    for (v in 1:V) {
      full_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[folds != v], X = draw$x[folds != v, , drop = FALSE], family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik"))
      full_fit[[v]] <- SuperLearner::predict.SuperLearner(full_mod, newdata = draw$x[folds == v, , drop = FALSE])$pred
      
      redu_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat[folds != v], X = draw$red_x[folds != v, , drop = FALSE], family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik")) # remove "method.NNloglik" due to numerical errors
      redu_fit[[v]] <- SuperLearner::predict.SuperLearner(redu_mod, newdata = draw$red_x[folds == v, , drop = FALSE])$pred
      
      full <- cbind(1 - full_fit[[v]], full_fit[[v]])
      redu <- cbind(1 - redu_fit[[v]], redu_fit[[v]])
      numerators[v] <- 2*sum(diag(t(full)%*%log(full/redu)))/dim(draw$y[folds == v, ])[1]
    }
    numerator <- mean(numerators)
  } else {
    full_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat, X = draw$x, family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik"))
    full_fit <- full_mod$SL.predict
    full <- cbind(1 - full_fit, full_fit) # predicting Y = 1
    red_mod <- suppressWarnings(SuperLearner::SuperLearner(Y = draw$y_cat, X = draw$red_x, family = binomial(), cvControl = list(V = V), SL.library = learner_lib, method = "method.CC_nloglik")) # remove "method.NNloglik" due to numerical errors
    red_fit <- red_mod$SL.predict
    redu <- cbind(1 - red_fit, red_fit)
    numerator <- 2*sum(diag(t(full)%*%log(full/redu)))/dim(draw$y)[1]
  }
  p <- apply(draw$y, 2, mean)
  denominator <- (-1)*sum(log(p))
  est <- numerator/denominator
  return(est)
}

## for bootstrapping
nonparametric_naive_deviance_boot <- function(data, indices, ...) {
  y <- data$y[indices, c(1, 2)]
  y_cat <- data$y_cat[indices]
  x <- data[indices, c(3, 4)]
  red_x <- data$red_x[indices]
  draw <- list(y = y, y_cat = y_cat, x = x, red_x = red_x)
  return(nonparametric_naive_deviance(draw, ...))
}
