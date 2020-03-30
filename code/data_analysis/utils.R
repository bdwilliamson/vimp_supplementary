#!/usr/bin/env Rscript

## plot vimp for a single outcome/group combination
## @param vimp_obj the variable importance object
## @param title the title of the plot
## @param x_lim the x-axis limits
## @param x_lab the x-axis label
## @param lgnd_pos the legend position
## @param point_size the point size
## @param main_font_size the size of text
## @param cv whether or not this is cv importance
## @param num_plot the number of features to plot (in descending order)
plot_one_vimp <- function(vimp_obj, title = "Variable importance", x_lim = c(0, 1), x_lab = expression(paste(R^2)), lgnd_pos = c(0.1, 0.3), point_size = 5, main_font_size = 20, plot_font_size = 6, axis_font_size = 18, cv = FALSE, num_plot = 50, threshold = 0.05, plot_vimp_nm = NULL) {
    text_pos <- x_lim[2] - 0.05
    signif_pos <- x_lim[2] - 0.025
    if (!is.null(vimp_obj)) {
        ## get the variable importances
        if (!is.null(vimp_obj$mat)) {
            vimp_est <- vimp_obj$mat
            vimp_est$group <- vimp_nice_rownames(vimp_obj, num_only = TRUE, cv = cv)
        } else {
            vimp_est <- cbind(est = vimp_obj$est, se = vimp_obj$se, cil = vimp_obj$cil, ciu = vimp_obj$ciu, p_value = vimp_obj$p_value)
            print_s <- ifelse(length(vimp_obj$s) <= 10,
                          paste(vimp_obj$s, collapse = ", "),
                          paste(c(vimp_obj$s[1:10], "..."), collapse = ", "))
            vimp_est$group <- paste("s = ", print_s, sep = "")
        }
        tmp_cis <- round(vimp_est[, c("cil", "ciu")], 3)
        text_cis <- apply(tmp_cis, 1, function(x) paste0("[", x[1], ", ", x[2], "]"))
        vimp_est <- tibble::add_column(vimp_est, text_ci = text_cis)
        vimp_est <- tibble::add_column(vimp_est, signif_p = ifelse(vimp_est$p_value < threshold, "*", ""))
        dim_plot <- ifelse(num_plot > dim(vimp_est)[1], dim(vimp_est)[1], num_plot)
        ## plot by ordered vimp measure
        vimp_plot <- vimp_est %>%
            arrange(desc(est)) %>%
            mutate(ord_group = forcats::fct_reorder(group, est)) %>%
            filter(row_number() <= num_plot) %>%
            ggplot(aes(x = est, y = ord_group)) +
            geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
            geom_point(size = point_size) +
            geom_text(aes(x = rep(text_pos, dim_plot), label = text_ci), hjust = "inward", size = plot_font_size, color = "blue") +
            geom_text(aes(x = rep(signif_pos, dim_plot), label = signif_p), hjust = "inward", size = plot_font_size, color = "blue") +
            ggtitle(title) +
            xlab(x_lab) +
            ylab("Feature group") +
            scale_x_continuous(breaks = round(seq(x_lim[1], x_lim[2], 0.1), 1),
                           labels = as.character(round(seq(x_lim[1], x_lim[2], 0.1), 1)),
                           limits = x_lim) +
            theme(axis.line = element_blank(), panel.border = element_rect(fill = NA, color = "black", linetype = 1, size = 1), plot.title = element_text(size = main_font_size), axis.title = element_text(size = axis_font_size), axis.text = element_text(size = axis_font_size), legend.text = element_text(size = axis_font_size), legend.title = element_text(size = axis_font_size), legend.position = lgnd_pos, text = element_text(size = main_font_size))
        return(vimp_plot)
    } else {
        print("No variable importance estimates to plot.")
    }
}

vimp_plot_name <- function(vimp_str) {
    plot_nms <- rep(NA, length(vimp_str))
    plot_nms[grepl("iip", vimp_str)] <- "IIP"
    plot_nms[grepl("pc.ic50", vimp_str)] <- "IC-50"
    plot_nms[grepl("pc.ic80", vimp_str)] <- "IC-80"
    plot_nms[grepl("dichotomous.1", vimp_str)] <- "Estimated sensitivity"
    plot_nms[grepl("dichotomous.2", vimp_str)] <- "Multiple sensitivity"
    plot_nms[grepl("censored", vimp_str)] <- "IC-50 Censored"
    plot_nms[grepl("sens.resis", vimp_str)] <- "Sensitive/Resistant"
    return(plot_nms)
}
vimp_plot_name_expr <- function(vimp_str) {
    plot_nms <- rep(NA, length(vimp_str))
    plot_nms[grepl("iip", vimp_str)] <- expression(bold(paste("IIP;")))
    plot_nms[grepl("pc.ic50", vimp_str)] <- expression(bold(paste(IC[50], ";")))
    plot_nms[grepl("pc.ic80", vimp_str)] <- expression(bold(paste(IC[80], ";")))
    plot_nms[grepl("dichotomous.1", vimp_str)] <- expression(bold(paste("Estimated sensitivity;")))
    plot_nms[grepl("dichotomous.2", vimp_str)] <- expression(bold(paste("Multiple sensitivity;")))
    plot_nms[grepl("censored", vimp_str)] <- expression(bold(paste(IC[50], " Censored;")))
    plot_nms[grepl("sens.resis", vimp_str)] <- expression(bold(paste("Sensitive/Resistant Only;")))
    return(plot_nms)
}
vimp_plot_name_expr_simple <- function(vimp_str) {
  plot_nms <- rep(NA, length(vimp_str))
  plot_nms[grepl("iip", vimp_str)] <- expression(bold(paste("IIP;")))
  plot_nms[grepl("pc.ic50", vimp_str)] <- expression(bold(paste(IC[50], ";")))
  plot_nms[grepl("pc.ic80", vimp_str)] <- expression(bold(paste(IC[80], ";")))
  plot_nms[grepl("dichotomous.1", vimp_str)] <- expression(bold(paste("Estimated sensitivity;")))
  plot_nms[grepl("dichotomous.2", vimp_str)] <- expression(bold(paste("Multiple sensitivity;")))
  plot_nms[grepl("censored", vimp_str)] <- expression(bold(paste("Resistant;")))
  plot_nms[grepl("sens.resis", vimp_str)] <- expression(bold(paste("Sensitive/Resistant Only;")))
  return(plot_nms)
}
vimp_nice_rownames <- function(vimp_obj, num_only = FALSE, cv = FALSE) {
    mat_s <- vimp_obj$mat$s
    lst_s <- vimp_obj$s
    indx_mat <- sapply(1:length(mat_s), function(x) which(mat_s[x] == lst_s))
    paste_ind <- 3
    if (cv) {
        paste_ind <- 4
    }
    tmp_nms <- unlist(lapply(strsplit(names(lst_s), "_", fixed = TRUE), function(x) paste(x[paste_ind:length(x)], collapse = "_")))
    if (num_only) {
        return(as.character(indx_mat))
    } else {
        return(tmp_nms[indx_mat])
    }
}

## ----------------------------------------------------
## SL functions
## ----------------------------------------------------
## plotting and summary functions for super learners with different methods

summary.myCV.SuperLearner <- function (object, obsWeights = NULL,
                                       method = NULL, ...) {
    if ("env" %in% names(object)) {
        env = object$env
    }
    else {
        env = parent.frame()
    }
    if(is.null(method)){
      method <- if (is.null(as.list(object$call)[["method"]])) {
        method <- "method.NNLS"
      }
      else if (is.symbol(as.list(object$call)[["method"]])) {
          method <- get(paste(as.list(object$call)[["method"]]),
              envir = env)
      }
      else {
          method <- as.list(object$call)[["method"]]
      }
    }
    library.names <- object$libraryNames
    V <- object$V
    n <- length(object$SL.predict)
    if (is.null(obsWeights)) {
        obsWeights <- rep(1, length(object$Y))
    }
    folds <- object$folds
    SL.predict <- object$SL.predict
    discreteSL.predict <- object$discreteSL.predict
    library.predict <- object$library.predict
    Y <- object$Y
    Risk.SL <- rep(NA, length = V)
    se.SL <- rep(NA, length = V)
    Risk.dSL <- rep(NA, length = V)
    se.dSL <- rep(NA, length = V)
    Risk.library <- matrix(NA, nrow = length(library.names),
        ncol = V)
    se.library <- matrix(NA, nrow = length(library.names),
        ncol = V)
    rownames(Risk.library) <- library.names
    if (method %in% c("method.NNLS", "method.NNLS2", "method.CC_LS")) {
        for (ii in seq_len(V)) {
            Risk.SL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                SL.predict[folds[[ii]]])^2)
            Risk.dSL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                discreteSL.predict[folds[[ii]]])^2)
            Risk.library[, ii] <- apply(library.predict[folds[[ii]],
                , drop = FALSE], 2, function(x) mean(obsWeights[folds[[ii]]] *
                (Y[folds[[ii]]] - x)^2))
        }
        if_sl <- (Y - SL.predict)^2 - mean((Y - SL.predict)^2)
        if_dsl <- (Y - discreteSL.predict)^2 - mean((Y - discreteSL.predict)^2)
        if_library <- apply(library.predict, 2, function(x){ (Y - x)^2 - mean((Y - x)^2) })
        if_varY <- (Y - mean(Y))^2 - mean((Y - mean(Y))^2)
        get_log_se <- function(if_risk, if_varY, risk, varY,
                               n = length(if_risk)){
            grad <- matrix(c(1 / risk, - 1 /varY), nrow = 2)
            Sig <- cov(cbind(if_risk, if_varY))
            se_log <- t(grad) %*% Sig %*% grad
            return(se_log)
        }

        se <- (1/sqrt(n)) * c(
          get_log_se(if_risk = if_sl, if_varY = if_varY, risk = mean(Risk.SL), varY = var(Y)),
          get_log_se(if_risk = if_dsl, if_varY = if_varY, risk = mean(Risk.dSL), varY = var(Y)),
          mapply(if1 = split(if_library, col(if_library)), risk = split(Risk.library, row(Risk.library)),
                 function(if1, risk){ get_log_se(if_risk = if1, if_varY = if_varY, risk = mean(risk), varY = var(Y))})
        )
    }
    else if (method %in% c("method.NNloglik", "method.CC_nloglik")) {
        for (ii in seq_len(V)) {
            Risk.SL[ii] <- -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]],
                log(SL.predict[folds[[ii]]]), log(1 - SL.predict[folds[[ii]]])))
            Risk.dSL[ii] <- -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]],
                log(discreteSL.predict[folds[[ii]]]), log(1 -
                  discreteSL.predict[folds[[ii]]])))
            Risk.library[, ii] <- apply(library.predict[folds[[ii]],
                , drop = FALSE], 2, function(x) {
                -mean(obsWeights[folds[[ii]]] * ifelse(Y[folds[[ii]]],
                  log(x), log(1 - x)))
            })
        }
        se <- rep.int(NA, (length(library.names) + 2))
    }
    else if (method %in% c("method.AUC")) {
        requireNamespace("cvAUC")
        for (ii in seq_len(V)) {
            sl_auc <- cvAUC::ci.cvAUC(predictions = SL.predict[folds[[ii]]],
                labels = Y[folds[[ii]]], folds = NULL)
            Risk.SL[ii] <- sl_auc$cvAUC
            se.SL[ii] <- sl_auc$se
            dsl_auc <- cvAUC::ci.cvAUC(predictions = discreteSL.predict[folds[[ii]]],
                labels = Y[folds[[ii]]], folds = NULL)
            Risk.dSL[ii] <- dsl_auc$cvAUC
            se.dSL[ii] <- dsl_auc$se
            library_auc <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x){
                tmp <- cvAUC::ci.cvAUC(predictions = x, labels = Y[folds[[ii]]], folds = NULL)
                return(c(tmp$cvAUC, tmp$se))
              })
            Risk.library[,ii] <- library_auc[1,]
            se.library[,ii] <- library_auc[2,]
        }
        se <- c(mean(se.SL, na.rm = TRUE), mean(se.dSL, na.rm = TRUE),
                rowMeans(se.library, na.rm = TRUE))
    }
    else {
        stop("summary function not available for SuperLearner with loss function/method used")
    }
    if(method != "method.AUC"){
      Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
          library.names), Ave = c(1 - mean(Risk.SL)/var(Y), 1 - mean(Risk.dSL)/var(Y),
          apply(Risk.library, 1, function(x){ 1 - mean(x)/var(Y) })), log_se = se, Min = c(min(1 - Risk.SL/var(Y)),
          min(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ min(1 - mean(x)/var(Y))})), Max = c(max(1 - Risk.SL/var(Y)),
          max(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ max(1 - mean(x)/var(Y)) })))
    }else{
      Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
        library.names), Ave = c(mean(Risk.SL), mean(Risk.dSL),
        apply(Risk.library, 1, mean)), se = se, Min = c(min(Risk.SL),
        min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(max(Risk.SL),
        max(Risk.dSL), apply(Risk.library, 1, max)))
    }
    out <- list(call = object$call, method = method, V = V, Risk.SL = Risk.SL,
          Risk.dSL = Risk.dSL, Risk.library = Risk.library, Table = Table)
    class(out) <- "summary.myCV.SuperLearner"
    return(out)
}

plot.myCV.SuperLearner <- function (x, package = "ggplot2", constant = qnorm(0.975), sort = TRUE, main_title = NULL,
                                    xlim1 = -0.025, xlim2 = 0.3, text_size = 10, Rsquared = TRUE,
    ...) {
    sumx <- summary(x, ...)
    if (sort)
        sumx$Table$Algorithm <- reorder(sumx$Table$Algorithm,
            sumx$Table$Ave)
    Mean <- sumx$Table$Ave
    if(Rsquared){
      se <- sumx$Table$log_se
      Lower <- 1 - exp( log(-Mean + 1) - constant * se)
      Upper <- 1 - exp( log(-Mean + 1) + constant * se)
    }else{
      se <- sumx$Table$se
      # put AUC CI on logit scale
      grad <- 1 / (Mean - Mean^2)
      logit_se <- sqrt(se^2 * grad^2)
      Lower <- plogis(qlogis(Mean) - constant * logit_se); Upper <- plogis(qlogis(Mean) + constant * logit_se)
    }
    assign("d", data.frame(Y = Mean, X = sumx$Table$Algorithm,
        Lower = Lower, Upper = Upper))
    if (package == "lattice") {
        .SL.require("lattice")
        p <- lattice::dotplot(X ~ Y, data = d, xlim = c(min(d$Lower) -
            0.02, max(d$Upper) + 0.02), xlab = "V-fold CV Risk Estimate",
            ylab = "Method", panel = function(x, y) {
                lattice::panel.xyplot(x, y, pch = 16, cex = 1)
                lattice::panel.segments(d$Lower, y, d$Upper,
                  y, lty = 1)
            })
    }
    if (package == "ggplot2") {
      this_y_lab <- if(Rsquared){
        expression("CV-"*R^2*" (95% CI)")
      }else{
        "CV-AUC"
      }
        SuperLearner:::.SL.require("ggplot2")
        p <- ggplot2::ggplot(d, ggplot2::aes_string(x = "X",
            y = "Y", ymin = "Lower", ymax = "Upper")) +
            ggplot2::geom_pointrange(size = 0.5) +
            ggplot2::coord_flip() + ggplot2::ylab(this_y_lab) +
            ggplot2::xlab("Method") +
            ggplot2::ggtitle(main_title) +
            scale_y_continuous(limits=c(xlim1, xlim2), oob = scales::squish) +
            ggplot2::theme(axis.text = element_text(size = text_size))
    }
    return(p)
}

plot_roc_curves <- function(cv_fit, topRank = 1,
                            cols = c(rgb(78, 103, 102, alpha = 255/2, maxColorValue = 255),
                                     rgb(90,177,187, alpha = 255/2, maxColorValue = 255),
                                     rgb(165, 200, 130, alpha = 255/2, maxColorValue = 255))) {
  class(cv_fit) <- "myCV.SuperLearner"
  allAlgos <- summary(cv_fit, method = "method.AUC")$Table %>% mutate(Algorithm = as.character(Algorithm))
  allCandidates <- allAlgos[-(1:2), ] %>% arrange(-Ave)
  sortedAlgos <- rbind(allAlgos[1:2,], allCandidates[1:topRank,])

  predict <- cv_fit[["library.predict"]] %>% as.data.frame() %>%
    bind_cols(cv_fit[["discreteSL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Discrete SL"))) %>%
    bind_cols(cv_fit[["SL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("Super Learner"))) %>%
    bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y"))) %>%
    gather("algo", "pred", -Y) %>%
    filter(algo %in% sortedAlgos$Algorithm[1:(2+topRank)])

  roc.obj <- predict %>%
    group_by(algo) %>%
    nest() %>%
    mutate(pred.obj = purrr::map(data, ~ ROCR::prediction(.x$pred, .x$Y)),
           perf.obj = purrr::map(pred.obj, ~ ROCR::performance(.x, "tpr", "fpr")),
           roc.dat = purrr::map(perf.obj, ~ tibble(xval = .x@x.values[[1]],
                                                   yval = .x@y.values[[1]])))
  roc.obj %>%
    unnest(roc.dat) %>%
    ggplot(aes(x=xval, y=yval, col=algo)) +
    geom_step(lwd=2) +
    theme(legend.position = "top") +
    scale_color_manual(values = cols) +
    labs(x = "Cross-Validated False Positive Rate", y = "Cross-Validated True Positive Rate", col = "Algorithm")
}
