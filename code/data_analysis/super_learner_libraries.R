#!/usr/local/bin/Rscript

## set up super learner library and optimization methods

## -------------------------------------------------------------------------
## candidate learners
## -------------------------------------------------------------------------
## boosted trees
SL.xgboost.corrected <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), nthread = 1, verbose = 0, save_period = NULL, ...) {
    SuperLearner:::.SL.require("xgboost")
    if (packageVersion("xgboost") < 0.6)
        stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
    if (!is.matrix(X)) {
        X = model.matrix(~. - 1, X)
    }
    xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
    if (family$family == "gaussian") {
        model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "binomial") {
        model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, nthread = nthread,
            params = params, save_period = save_period)
    }
    if (family$family == "multinomial") {
        model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
            nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
            eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
            nthread = nthread, params = params, save_period = save_period)
    }
    if (!is.matrix(newX)) {
        newX = model.matrix(~. - 1, newX)
    }
    pred = predict(model, newdata = newX)
    fit = list(object = model)
    class(fit) = c("SL.xgboost")
    out = list(pred = pred, fit = fit)
    return(out)
}
## xgboost with different max_depth parameters
SL.xgboost1 <- function(..., max_depth = 1) {
    SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost2 <- function(..., max_depth = 2){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost4 <- function(..., max_depth = 4){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost6 <- function(..., max_depth = 6){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost8 <- function(..., max_depth = 8){
	SL.xgboost.corrected(..., max_depth = max_depth)
}
descr_SL.xgboost <- "boosted regression trees with maximum depth of "
descr_SL.xgboost1 <- paste0(descr_SL.xgboost, 1)
descr_SL.xgboost2 <- paste0(descr_SL.xgboost, 2)
descr_SL.xgboost4 <- paste0(descr_SL.xgboost, 4)
descr_SL.xgboost6 <- paste0(descr_SL.xgboost, 6)
descr_SL.xgboost8 <- paste0(descr_SL.xgboost, 8)
# random forests
SL.ranger.imp <- function (Y, X, newX, family, obsWeights, num.trees = 500, mtry = floor(sqrt(ncol(X))), write.forest = TRUE, probability = family$family == "binomial", min.node.size = ifelse(family$family == "gaussian", 5, 1), replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632), num.threads = 1, verbose = TRUE, ...) {
    SuperLearner:::.SL.require("ranger")
    if (family$family == "binomial") {
        Y = as.factor(Y)
    }
    if (is.matrix(X)) {
        X = data.frame(X)
    }
    fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
        replace = replace, sample.fraction = sample.fraction,
        case.weights = obsWeights, write.forest = write.forest,
        probability = probability, num.threads = num.threads,
        verbose = verbose, importance = "impurity")
    pred <- predict(fit, data = newX)$predictions
    if (family$family == "binomial") {
        pred = pred[, "1"]
    }
    fit <- list(object = fit, verbose = verbose)
    class(fit) <- c("SL.ranger")
    out <- list(pred = pred, fit = fit)
    return(out)
}
SL.ranger.reg <- function(..., X, mtry = floor(sqrt(ncol(X)))){
	SL.ranger.imp(..., X = X, mtry = mtry)
}

SL.ranger.small <- function(..., X, mtry = floor(sqrt(ncol(X)) * 1/2)){
	SL.ranger.imp(..., X = X, mtry = mtry)
}

SL.ranger.large <- function(..., X, mtry = floor(sqrt(ncol(X)) * 2)){
	SL.ranger.imp(..., X = X, mtry = mtry)
}
descr_SL.ranger.imp <- "random forest with mtry equal to "
descr_SL.ranger.reg <- paste0(descr_SL.ranger.imp, "square root of number of predictors")
descr_SL.ranger.small <- paste0(descr_SL.ranger.imp, "one-half times square root of number of predictors")
descr_SL.ranger.large <- paste0(descr_SL.ranger.imp, "two times square root of number of predictors")

# function used to do smarter CV for glmnet
get_fold_id <- function(Y){
  fold_id <- rep(0, length(Y))
  wiY0 <- which(Y == 0)
  wiY1 <- which(Y == 1)
  #if <4 cases, no cv
  if(length(wiY1) == 4){
    #if exactly 4 cases, 4-fold cv
    #1 case per fold
    fold <- 1:4
    fold_id[sample(wiY1)] <- fold
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }else{
    #if >=5 cases, 5 fold cv
    #cases split as evenly as possible
    fold <- 1:5
    fold_id[sample(wiY1)] <- rep(fold, length = length(wiY1))
    fold_id[sample(wiY0)] <- rep(fold, length = length(wiY0))
  }
  return(fold_id)
}


# function to have more robust behavior in SL.glmnet
SL.glmnet.mycv <- function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 5,
    nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
    SuperLearner:::.SL.require("glmnet")
    if (!is.matrix(X)) {
        X <- model.matrix(~-1 + ., X)
        newX <- model.matrix(~-1 + ., newX)
    }
    fold_id <- get_fold_id(Y)
    nfolds <- max(fold_id)
    if(nfolds != 0){
        fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
            lambda = NULL, type.measure = loss, nfolds = nfolds,
            family = family$family, alpha = alpha, nlambda = nlambda,
            ...)
        pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
            "lambda.min", "lambda.1se"))
        fit <- list(object = fitCV, useMin = useMin)
        class(fit) <- "SL.glmnet"
    }else{
        # if fewer than 3 cases, just use mean
        meanY <- weighted.mean(Y, w = obsWeights)
        pred <- rep.int(meanY, times = nrow(newX))
        fit <- list(object = meanY)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.mean")
    }
    out <- list(pred = pred, fit = fit)
    return(out)
}

# elastic net with different levels of alpha
SL.glmnet.1 <- function(..., alpha = 1) {
    SL.glmnet.mycv(..., alpha = alpha)
}
SL.glmnet.50 <- function(..., alpha = 0.5){
	SL.glmnet.mycv(..., alpha = alpha)
}
SL.glmnet.25 <- function(..., alpha = 0.25){
	SL.glmnet.mycv(..., alpha = alpha)
}
SL.glmnet.75 <- function(..., alpha = 0.75){
	SL.glmnet.mycv(..., alpha = alpha)
}

descr_SL.glmnet <- "GLMNET with lambda selected by 5-fold CV and alpha equal to "
descr_SL.glmnet.50 <- paste0(descr_SL.glmnet, "0.5")
descr_SL.glmnet.25 <- paste0(descr_SL.glmnet, "0.25")
descr_SL.glmnet.75 <- paste0(descr_SL.glmnet, "0.75")
descr_SL.glmnet.1 <- paste0(descr_SL.glmnet, "1")
descr_SL.glmnet.mycv <- "GLMNET with lambda selected by CV and alpha equal to 0"

descr_SL.mean <- "intercept only regression"
descr_SL.glm <- "main terms generalized linear model"

library_full <- c("SL.mean", paste0("SL.xgboost", c(1, 2, 4, 6, 8)), paste0("SL.ranger.", c("small", "reg", "large")), "SL.glmnet.mycv", paste0("SL.glmnet.", c(25, 50, 75, 1)))

library_reduced <- c("SL.ranger.reg", "SL.glmnet.1", "SL.xgboost1")

## -------------------------------------------------------------------------
## temporary optimization methods
## -------------------------------------------------------------------------
#' Temporary fix for convex combination method mean squared error
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another
tmp_method.CC_LS <- function () {
    computeCoef = function(Z, Y, libraryNames, verbose, obsWeights,
        errorsInLibrary = NULL, ...) {
        cvRisk <- apply(Z, 2, function(x) mean(obsWeights * (x -
            Y)^2))
        names(cvRisk) <- libraryNames
        compute <- function(x, y, wt = rep(1, length(y))) {
            wX <- sqrt(wt) * x
            wY <- sqrt(wt) * y
            D <- crossprod(wX)
            d <- crossprod(wX, wY)
            A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
            bvec <- c(1, rep(0, ncol(wX)))
            fit <- tryCatch({quadprog::solve.QP(Dmat = D, dvec = d, Amat = A,
                bvec = bvec, meq = 1)
          }, error = function(e){
            out <- list()
            class(out) <- "error"
            out
          })
            invisible(fit)
        }
        modZ <- Z
        naCols <- which(apply(Z, 2, function(z) {
            all(z == 0)
        }))
        anyNACols <- length(naCols) > 0
        if (anyNACols) {
            warning(paste0(paste0(libraryNames[naCols], collapse = ", "),
                " have NAs.", "Removing from super learner."))
        }
        tol <- 4
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "),
                " are duplicates of previous learners.", " Removing from super learner."))
        }
        if (anyDupCols | anyNACols) {
            rmCols <- unique(c(naCols, dupCols))
            modZ <- Z[, -rmCols, drop = FALSE]
        }
        fit <- compute(x = modZ, y = Y, wt = obsWeights)
        if(class(fit) != "error"){
          coef <- fit$solution
        }else{
          coef <- rep(0, ncol(Z))
          coef[which.min(cvRisk)] <- 1
        }
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] = 0
        }
        if(class(fit) != "error"){
          if (anyDupCols | anyNACols) {
              ind <- c(seq_along(coef), rmCols - 0.5)
              coef <- c(coef, rep(0, length(rmCols)))
              coef <- coef[order(ind)]
          }
          coef[coef < 1e-04] <- 0
          coef <- coef/sum(coef)
        }
        if (!sum(coef) > 0)
            warning("All algorithms have zero weight", call. = FALSE)
        list(cvRisk = cvRisk, coef = coef, optimizer = fit)
    }
    computePred = function(predY, coef, ...) {
        predY %*% matrix(coef)
    }
    out <- list(require = "quadprog", computeCoef = computeCoef,
        computePred = computePred)
    invisible(out)
}


#' Temporary fix for convex combination method negative log-likelihood loss
#' Relative to existing implementation, we reduce the tolerance at which
#' we declare predictions from a given algorithm the same as another.
#' Note that because of the way \code{SuperLearner} is structure, one needs to
#' install the optimization software separately.
tmp_method.CC_nloglik <- function () {
    computePred = function(predY, coef, control, ...) {
        if (sum(coef != 0) == 0) {
            stop("All metalearner coefficients are zero, cannot compute prediction.")
        }
        stats::plogis(trimLogit(predY[, coef != 0], trim = control$trimLogit) %*%
            matrix(coef[coef != 0]))
    }
    computeCoef = function(Z, Y, libraryNames, obsWeights, control,
        verbose, ...) {
        tol <- 4
        dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
        anyDupCols <- length(dupCols) > 0
        modZ <- Z
        if (anyDupCols) {
            warning(paste0(paste0(libraryNames[dupCols], collapse = ", "),
                " are duplicates of previous learners.", " Removing from super learner."))
            modZ <- modZ[, -dupCols, drop = FALSE]
        }
        modlogitZ <- trimLogit(modZ, control$trimLogit)
        logitZ <- trimLogit(Z, control$trimLogit)
        cvRisk <- apply(logitZ, 2, function(x) -sum(2 * obsWeights *
            ifelse(Y, stats::plogis(x, log.p = TRUE), stats::plogis(x, log.p = TRUE,
                lower.tail = FALSE))))
        names(cvRisk) <- libraryNames
        obj_and_grad <- function(y, x, w = NULL) {
            y <- y
            x <- x
            function(beta) {
                xB <- x %*% cbind(beta)
                loglik <- y * stats::plogis(xB, log.p = TRUE) + (1 -
                  y) * stats::plogis(xB, log.p = TRUE, lower.tail = FALSE)
                if (!is.null(w))
                  loglik <- loglik * w
                obj <- -2 * sum(loglik)
                p <- stats::plogis(xB)
                grad <- if (is.null(w))
                  2 * crossprod(x, cbind(p - y))
                else 2 * crossprod(x, w * cbind(p - y))
                list(objective = obj, gradient = grad)
            }
        }
        lower_bounds = rep(0, ncol(modZ))
        upper_bounds = rep(1, ncol(modZ))
        if (anyNA(cvRisk)) {
            upper_bounds[is.na(cvRisk)] = 0
        }
        r <- tryCatch({nloptr::nloptr(x0 = rep(1/ncol(modZ), ncol(modZ)),
            eval_f = obj_and_grad(Y, modlogitZ), lb = lower_bounds,
            ub = upper_bounds, eval_g_eq = function(beta) (sum(beta) -
                1), eval_jac_g_eq = function(beta) rep(1, length(beta)),
            opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_abs = 1e-08))
        }, error = function(e){
          out <- list()
          class(out) <- "error"
          out
        })
        if (r$status < 1 || r$status > 4) {
            warning(r$message)
        }
        if(class(r) != "error"){
          coef <- r$solution
        }else{
          coef <- rep(0, ncol(Z))
          coef[which.min(cvRisk)] <- 1
        }
        if (anyNA(coef)) {
            warning("Some algorithms have weights of NA, setting to 0.")
            coef[is.na(coef)] <- 0
        }
        if (anyDupCols) {
            ind <- c(seq_along(coef), dupCols - 0.5)
            coef <- c(coef, rep(0, length(dupCols)))
            coef <- coef[order(ind)]
        }
        coef[coef < 1e-04] <- 0
        coef <- coef/sum(coef)
        out <- list(cvRisk = cvRisk, coef = coef, optimizer = r)
        return(out)
    }
    list(require = "nloptr", computeCoef = computeCoef, computePred = computePred)
}
