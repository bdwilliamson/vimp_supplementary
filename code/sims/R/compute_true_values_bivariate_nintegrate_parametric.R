## compute true values by numerical integration for simple data-generating mechanism

# make binomial y's
make_y <- function(b, p) rbinom(b, 1, p)
## make independent x's
make_x_parametric <- function(b, mu_0, mu_1, sigma, y) {
  n_col <- ifelse(!is.null(dim(mu_0)), dim(mu_0)[2], length(mu_0))
  x <- matrix(0, nrow = b, ncol = n_col)
  x[y == 0, ] <- MASS::mvrnorm(n = sum(y == 0), mu = mu_0, Sigma = sigma)
  x[y == 1, ] <- MASS::mvrnorm(n = sum(y == 1), mu = mu_1, Sigma = sigma)
  if (n_col == 1) {
    x <- cbind(x, rnorm(b, 0, 1))
  }
  return(x)
}

# put a 1 in the correct column of a matrix, given a number
make_mat <- function(x, K) {
  vec <- matrix(0, nrow = 1, ncol = K)
  vec[x+1] <- 1
  return(vec)
}
# turn a k-class random variable into a k-class multinomial random vector
create_multinomial <- function(x) {
  # get the number of classes
  K <- length(unique(x))
  # create a length K vector for each observation
  # with a 1 in position k
  ret_mat <- apply(matrix(x), 1, function(x) make_mat(x, K))
}

gen_data_parametric <- function(a, p, mu_0, mu_1, sigma) {
  # create y
  y <- make_y(a, p)
  # create x
  x <- make_x_parametric(a, mu_0, mu_1, sigma, y)
  # make categorical
  y_cat <- y
  y <- t(create_multinomial(y_cat))
  return(list(x = x, y = y, y_cat = y_cat))
}
p <- 0.6
n <- 1000
mu_0 <- matrix(c(0, 0), nrow = 1)
mu_1 <- matrix(c(1.5, 2), nrow = 1)
mu_0_mat <- matrix(rep(mu_0, n), nrow = n, byrow = TRUE)
mu_1_mat <- matrix(rep(mu_1, n), nrow = n, byrow = TRUE)
Sigma <- diag(1, nrow = 2)
# outer_nsim <- 10000
outer_nsim <- 100000

## full function of x
f <- function(x, mu_0, mu_1, Sigma, p_y) {
    numerator <- exp((-1/2)*(x - mu_1)%*%solve(Sigma)%*%t(x - mu_1))*p_y
    denominator <- numerator + exp((-1/2)*(x - mu_0)%*%solve(Sigma)%*%t(x - mu_0))*(1 - p_y)
    return(numerator/denominator)
}
## reduced function is simply with smaller matrices (remove the correct position)

## denominator for all deviance calculations
deviance_denom <- (-1)*(log(p) + log(1 - p))

## loop through to get deviance numerator, etc. for each sim
set.seed(4747)
deviance_nums_1 <- deviance_nums_2 <- adas_1 <- adas_2 <- daucs_1 <- daucs_2 <- vector("numeric", outer_nsim)
system.time(for (i in 1:outer_nsim) {
    ## generate data
    data <- gen_data_parametric(n, p, mu_0, mu_1, Sigma)
    ## compute full, reduced for the given dataset; for each obs (in rows)
    full <- diag(f(data$x, mu_0_mat, mu_1_mat, Sigma, p))
    reduced_remove_1 <- diag(f(data$x[, -1, drop = FALSE], mu_0_mat[, -1, drop = FALSE], mu_1_mat[, -1, drop = FALSE], diag(1), p))
    reduced_remove_2 <- diag(f(data$x[, -2, drop = FALSE], mu_0_mat[, -2, drop = FALSE], mu_1_mat[, -2, drop = FALSE], diag(1), p))
    full_mat <- cbind(1 - full, full)
    reduced_mat_1 <- cbind(1 - reduced_remove_1, reduced_remove_1)
    reduced_mat_2 <- cbind(1 - reduced_remove_2, reduced_remove_2)
    ## compute deviance
    deviance_nums_1[i] <- 2*sum(diag(t(data$y)%*%log(full_mat/reduced_mat_1)), na.rm = TRUE)/dim(data$y)[1]
    deviance_nums_2[i] <- 2*sum(diag(t(data$y)%*%log(full_mat/reduced_mat_2)), na.rm = TRUE)/dim(data$y)[1]

    ## compute ada
    contrib_ada_full <- mean((full > 1/2) != data$y_cat)
    contrib_ada_rem1 <- mean((reduced_remove_1 > 1/2) != data$y_cat)
    contrib_ada_rem2 <- mean((reduced_remove_2 > 1/2) != data$y_cat)
    adas_1[i] <- (1 - contrib_ada_full) - (1 - contrib_ada_rem1)
    adas_2[i] <- (1 - contrib_ada_full) - (1 - contrib_ada_rem2)

    ## compute dauc
    full_pred <- ROCR::prediction(predictions = full, labels = data$y_cat)
    red1_pred <- ROCR::prediction(predictions = reduced_remove_1, labels = data$y_cat)
    red2_pred <- ROCR::prediction(predictions = reduced_remove_2, labels = data$y_cat)

    auc_full <- unlist(ROCR::performance(prediction.obj = full_pred, measure = "auc", x.measure = "cutoff")@y.values)
    auc_rem1 <- unlist(ROCR::performance(prediction.obj = red1_pred, measure = "auc", x.measure = "cutoff")@y.values)
    auc_rem2 <- unlist(ROCR::performance(prediction.obj = red2_pred, measure = "auc", x.measure = "cutoff")@y.values)
    daucs_1[i] <- auc_full - auc_rem1
    daucs_2[i] <- auc_full - auc_rem2
})
deviance_1 <- mean(deviance_nums_1)/deviance_denom
deviance_2 <- mean(deviance_nums_2)/deviance_denom

ada_1 <- mean(adas_1)
ada_2 <- mean(adas_2)

dauc_1 <- mean(daucs_1)
dauc_2 <- mean(daucs_2)

truths <- data.frame(j_1 = c(deviance_1, ada_1, dauc_1), 
                     j_2 = c(deviance_2, ada_2, dauc_2), 
                     type = c("deviance", "accuracy", "auc"))

## save it 
saveRDS(truths, "truths_binary_bivariate.rds")

## loop through to get deviance numerator, etc. for each sim; under null
mu_2 <- 0
mu_0_null <- matrix(c(0, mu_2), nrow = 1)
mu_1_null <- matrix(c(1.5, mu_2), nrow = 1)
mu_0_null_mat <- matrix(rep(mu_0_null, n), nrow = n, byrow = TRUE)
mu_1_null_mat <- matrix(rep(mu_1_null, n), nrow = n, byrow = TRUE)
Sigma_null <- diag(1, nrow = 2)

deviance_null_nums_1 <- adas_null_1 <- daucs_null_1 <- vector("numeric", outer_nsim)
set.seed(5555)
system.time(for (i in 1:outer_nsim) {
  ## generate data
  data <- gen_data_parametric(n, p, mu_0_null, mu_1_null, Sigma_null)
  ## compute full, reduced for the given dataset; for each obs (in rows)
  full <- diag(f(data$x, mu_0_null_mat, mu_1_null_mat, Sigma_null, p))
  reduced_remove_1 <- diag(f(data$x[, -1, drop = FALSE], mu_0_null_mat[, -1, drop = FALSE], mu_1_null_mat[, -1, drop = FALSE], diag(1), p))
  # reduced_remove_1 <- rep(p, n)
  full_mat <- cbind(1 - full, full)
  reduced_mat_1 <- cbind(1 - reduced_remove_1, reduced_remove_1)
  ## compute deviance
  deviance_num <- 2*sum(diag(t(data$y)%*%log(full_mat/reduced_mat_1)), na.rm = TRUE)/dim(data$y)[1]
  deviance_null_nums_1[i] <- deviance_num

  ## compute ada
  contrib_ada_full <- mean((full > 1/2) != data$y_cat)
  contrib_ada_rem1 <- mean((reduced_remove_1 > 1/2) != data$y_cat)
  adas_null_1[i] <- (1 - contrib_ada_full) - (1 - contrib_ada_rem1)

  ## compute dauc
  full_pred <- ROCR::prediction(predictions = full, labels = data$y_cat)
  # redu_pred <- ROCR::prediction(predictions = reduced_remove_1, labels = data$y_cat)

  auc_full <- unlist(ROCR::performance(prediction.obj = full_pred, measure = "auc", x.measure = "cutoff")@y.values)
  # auc_redu <- unlist(ROCR::performance(prediction.obj = redu_pred, measure = "auc", x.measure = "cutoff")@y.values)
  daucs_null_1[i] <- auc_full - 0.5
})
deviance_null_1 <- mean(deviance_null_nums_1)/deviance_denom
deviance_null_2 <- 0

ada_null_1 <- mean(adas_null_1)
ada_null_2 <- 0

auc_null_1 <- mean(daucs_null_1) 
auc_null_2 <- 0

truths_null <- data.frame(j_1 = c(deviance_null_1, ada_null_1, auc_null_1), 
                          j_2 = c(deviance_null_2, ada_null_2, auc_null_2), 
                          type = c("deviance", "accuracy", "auc"))
saveRDS(truths_null, "truths_binary_bivariate_null.rds")
