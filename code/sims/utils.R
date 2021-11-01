# utility functions

# get the true values
# @param args the list of arguments
# @param the number of digits to round to
# @return the true values
get_true_values <- function(args, num_digits, dir) {
  if (args$p > 4) {
    noise_truths <- tibble::tibble(type = c("deviance", "accuracy", "auc")) %>%
      bind_cols(
        as.data.frame(matrix(0, nrow = 3,
                             ncol = args$p - 4,
                             dimnames = list(NULL,
                                             paste0("j_",
                                                    5:args$p))))
      )
  } else {
    noise_truths <- tibble::tibble(type = c("deviance", "accuracy", "auc"))
  }
  if (args$corr > 0) {
    truths <- readRDS(paste0(dir, "/true_vals_probit_corr.rds")) %>%
      left_join(noise_truths, by = "type")
  } else {
    truths <- readRDS(paste0(dir, "/true_vals_probit.rds")) %>%
      left_join(noise_truths, by = "type")
  }
  truths %>%
    select(starts_with("j_"), "type") %>%
    mutate(across(starts_with("j_"), round, digits = num_digits)) %>%
    filter(type != "deviance")
}

# get SL options
# @param y the outcome
# @param learners the learners
# @return the family, method, and learners for SL
get_sl_opts <- function(y = NULL, link_fun = "logit", learners = c("SL.gam", "SL.ranger")) {
  if (length(unique(y)) == 2) {
    family <- binomial(link = link_fun)
    method <- "method.CC_nloglik"
  } else {
    family <- gaussian()
    method <- "method.CC_LS"
  }
  list(fam = family, method = method, learners = learners)
}
# get ranger tuning parameters
# @param n the sample size
# @param family the type of regression to run
# @return the max depth, number of trees, min node size, probability, etc. for ranger
get_rf_opts <- function(draw, family = binomial(link = "probit"),
                        num_trees = c(100, 500), min_node_size = c(1, 3, 5),
                        max_depth = 1) {
  dat <- cbind.data.frame(Y = factor(draw$y_cat), as.data.frame(draw$x))
  folds <- vimp::make_folds(y = draw$y_cat, V = 5, stratified = TRUE, C = rep(1, nrow(draw$x)))
  grid <- expand.grid(ntree = num_trees, min_node = min_node_size, max_depth = max_depth,
                      v = 1:5)
  grid$log_loss <- vector("numeric", length = nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fit <- ranger::ranger(Y ~ ., data = subset(dat, folds != grid$v[i]),
                          replace = TRUE, probability = (family$family == "binomial"),
                          num.threads = 1, sample.fraction = 1,
                          num.trees = grid$ntree[i], min.node.size = grid$min_node[i],
                          max.depth = grid$max_depth[i],
                          mtry = switch(
                            ((ncol(dat) - 1) == 1) + 1, floor(sqrt(ncol(dat))),
                            ncol(dat) - 1
                          ))
    preds <- predict(fit, data = subset(dat, folds == grid$v[i]))$predictions[, "1"]
    test_y <- draw$y_cat[folds == grid$v[i]]
    grid$log_loss[i] <- (-1) * mean( test_y * log(preds) + (1 - test_y) * log(1 - preds) , na.rm = TRUE)
  }
  cv_perf <- grid %>%
    group_by(ntree, min_node, max_depth) %>%
    summarize(mn_log_loss = mean(log_loss), .groups = "drop")
  indx <- which.min(cv_perf$mn_log_loss)
  list(num_trees = cv_perf$ntree[indx], min_node_size = cv_perf$min_node[indx],
       max_depth = cv_perf$max_depth[indx])
}
# get xgboost tuning parameters
# @param n the sample size
# @param family the type of regression to run
# @return
get_xgb_opts <- function(draw, family = binomial(link = "probit"),
                         shrinkage = c(0.01, 0.1), max_depth = c(2, 4)) {
  dat <- cbind.data.frame(Y = draw$y_cat, as.data.frame(draw$x))
  folds <- vimp::make_folds(y = draw$y_cat, V = 5, stratified = TRUE, C = rep(1, nrow(draw$x)))
  objective <- switch((family$family == "binomial") + 1, "reg:squarederror",
                      "binary:logistic")
  grid <- expand.grid(shrinkage = shrinkage, max_depth = max_depth, v = 1:5)
  grid$log_loss <- vector("numeric", length = nrow(grid))
  for (i in seq_len(nrow(grid))) {
    train_df <- subset(dat, folds != grid$v[i])
    test_df <- subset(dat, folds == grid$v[i])
    train_mat <- model.matrix(~ . - 1, train_df[, -1, drop = FALSE])
    test_mat <- model.matrix(~ . - 1, test_df[, -1, drop = FALSE])
    xgb_train <- xgboost::xgb.DMatrix(data = train_mat, label = train_df$Y)
    fit <- xgboost::xgboost(data = xgb_train, objective = objective,
                            params = list(eta = grid$shrinkage[i],
                                          max_depth = grid$max_depth[i]),
                            nrounds = 1000, verbose = 0)
    preds <- predict(fit, newdata = test_mat)
    test_y <- draw$y_cat[folds == grid$v[i]]
    grid$log_loss[i] <- (-1) * mean( test_y * log(preds) + (1 - test_y) * log(1 - preds) , na.rm = TRUE)
  }
  cv_perf <- grid %>%
    group_by(shrinkage, max_depth) %>%
    summarize(mn_log_loss = mean(log_loss), .groups = "drop")
  indx <- which.min(cv_perf$mn_log_loss)
  list(eta = cv_perf$shrinkage[indx], max_depth = cv_perf$max_depth[indx])
}

# Get a list of reduced data objects for sample-splitting
# @param draw the original data
# @param ss_folds the sample-splitting folds (for the full data)
# @param full is this the full or reduced fit?
get_sample_split_data <- function(draw, ss_folds, full = TRUE) {
  new_draw <- list()
  split <- ifelse(full, 1, 2)
  new_draw$y <- draw$y[ss_folds == split, ]
  new_draw$y_cat <- draw$y_cat[ss_folds == split] 
  new_draw$x <- subset(draw$x, ss_folds == split) 
  new_draw$red_x <- subset(draw$red_x, ss_folds == split) 
  new_draw$j <- draw$j
  new_draw
}
