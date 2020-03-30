#!/usr/local/bin/Rscript
## ------------------------------------------------
## FILE: sim_binary_bivariate_data.R
## CREATED: 06 March 2019 by Brian Williamson
## PURPOSE: create one dataset for binary sim, 
##          bivariate predictor
## ------------------------------------------------

# expit of a value
expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}
logit <- function(x) {
  return(log(x/(1-x)))
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

# make binomial y's
make_y <- function(b, p) rbinom(b, 1, p)
## make independent x's
# make_x <- function(b) {
#   x <- replicate(2, rnorm(b, mean = 0, sd = 1))
#   return(x)
# }
make_x <- function(b, mu_0, mu_1, sigma, y) {
    n_col <- ifelse(!is.null(dim(mu_0)), dim(mu_0)[2], length(mu_0))
    x <- matrix(0, nrow = b, ncol = n_col)
    x[y == 0, ] <- MASS::mvrnorm(n = sum(y == 0), mu = mu_0, Sigma = sigma)
    x[y == 1, ] <- MASS::mvrnorm(n = sum(y == 1), mu = mu_1, Sigma = sigma)
    if (n_col == 1) {
      x <- cbind(x, rnorm(b, 0, 1))
    }
    return(x)
}
## generate a dataset
gen_data <- function(a, mu_0, mu_1, sigma, p, j) {
  ## create y
  y <- make_y(a, p)
  # create x
  # x <- make_x(a)
  x <- make_x(a, mu_0, mu_1, sigma, y)
  red_x <- x[, -j]
  
  # make categorical
  y_cat <- y
  y <- t(create_multinomial(y_cat))
  return(list(x = x, red_x = red_x, y = y, y_cat = y_cat, j = j))
}
