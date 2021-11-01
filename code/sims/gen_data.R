#!/usr/local/bin/Rscript
# generate data for simulations

# put a 1 in the correct column of a matrix, given a number
make_mat <- function(x = NULL, K = 2) {
  vec <- matrix(0, nrow = 1, ncol = K)
  vec[x+1] <- 1
  return(vec)
}
# turn a k-class random variable into a k-class multinomial random vector
create_multinomial <- function(x = NULL) {
  # get the number of classes
  K <- length(unique(x))
  # create a length K vector for each observation
  # with a 1 in position k
  ret_mat <- apply(matrix(x), 1, function(x) make_mat(x, K))
}

# make binomial y values
make_y <- function(b = 100, p = rep(0.5, b)) rbinom(b, 1, p)
# make x matrix:
#   if p = 2, then no correlation
#   if p > 2 and rho = 0, no correlation
#   if p > 2 and rho != 0, correlation
make_x <- function(b = 100, mu_0 = matrix(c(1, 1)), mu_1 = matrix(c(1, 1)), sigma = 0.25, y = rep(1, b), p = 2, rho = c(0, 0)) {
    n_col <- ifelse(!is.null(dim(mu_0)), dim(mu_0)[2], length(mu_0))
    x <- matrix(0, nrow = b, ncol = n_col)
    if (b == 1) {
      if (y == 0) {
        x <- matrix(MASS::mvrnorm(n = 1, mu = mu_0, Sigma = sigma), nrow = 1)
      } else {
        x <- matrix(MASS::mvrnorm(n = 1, mu = mu_1, Sigma = sigma), nrow = 1)
      }
    } else {
      x[y == 0, ] <- MASS::mvrnorm(n = sum(y == 0), mu = mu_0, Sigma = sigma)
      x[y == 1, ] <- MASS::mvrnorm(n = sum(y == 1), mu = mu_1, Sigma = sigma)
    }
    if (n_col == 1) {
      x <- cbind(x, rnorm(b, 0, 1))
    }
    if (p > 2 && ncol(x) < p) {
        if (any(rho == 0)) {
            noise_x <- t(replicate(p - 2, rnorm(b, 0, 1)))
        } else {
            noise_x <- t(replicate(p - 4, rnorm(b, 0, 1)))
        }
        x <- cbind(x, noise_x)
    }
    return(x)
}
# generate a dataset
gen_logit_data <- function(a = 100, mu_0 = matrix(c(1, 1)), mu_1 = matrix(c(1, 1)), 
                           sigma = 0.25, prob_y = 0.6, j = 1, p = 2, rho = 0) {
  # create y
  y <- make_y(b = a, p = prob_y)
  # create x
  x <- make_x(b = a, mu_0 = mu_0, mu_1 = mu_1, sigma = sigma, y = y, p = p, rho = rho)
  red_x <- x[, -j]

  # make categorical
  y_cat <- y
  y <- t(create_multinomial(y_cat))
  return(list(x = x, red_x = red_x, y = y, y_cat = y_cat, j = j))
}
# generate probit data
gen_probit_data <- function(nsim = 100, p = 4, x_mean = rep(0, p), Sigma = diag(1, p),
                     beta_0 = matrix(c(1.5, 2, rep(0, p - 2))), j = 1) {
  x <- MASS::mvrnorm(n = nsim, mu = x_mean, Sigma = Sigma)
  y <- as.numeric((x %*% beta_0 + rnorm(nsim)) > 0)
  red_x <- x[, -j, drop = FALSE]
  y_cat <- y
  y <- t(create_multinomial(y_cat))
  return(list(x = x, red_x = red_x, y = y, y_cat = y_cat, j = j))
}