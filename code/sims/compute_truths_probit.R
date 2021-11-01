# compute truths under a probit model
library("here")

# set up the data
sample_size <- 1e6
p <- 4
x_mean <- rep(0, p)
beta_0 <- matrix(c(2.5, 3.5, 0, 0))
Sigma <- diag(1, p)

# ------------------------------------------------------------------------------
# CORRELATED DATA
# ------------------------------------------------------------------------------
rho_13 <- 0.7
rho_24 <- 0.2
Sigma_corr <- Sigma
Sigma_corr[1, 3] <- Sigma_corr[3, 1] <- rho_13
Sigma_corr[2, 4] <- Sigma_corr[4, 2] <- rho_24

set.seed(4747)
x <- MASS::mvrnorm(n = sample_size, mu = x_mean, Sigma = Sigma_corr)
y <- as.numeric((x %*% beta_0 + rnorm(sample_size)) > 0)
# compute the various population prediction functions
f_0 <- pnorm(x %*% beta_0)
f_01 <- pnorm((x[, 2] * beta_0[2] + beta_0[1] * rho_13 * x[, 3]) /
                sqrt(1 + beta_0[1] ^ 2 * (1 - rho_13 ^ 2)))
f_02 <- pnorm((x[, 1] * beta_0[1] + beta_0[2] * rho_24 * x[, 4]) /
                sqrt(1 + beta_0[2] ^ 2 * (1 - rho_24 ^ 2)))
f_03 <- f_04 <- f_0
f_013 <- pnorm(beta_0[2] * x[, 2])
f_024 <- pnorm(beta_0[1] * x[, 1])

# compute the variable importance estimates
full_accuracy <- vimp::measure_accuracy(f_0, y)$point_est
acc_vim_1 <- full_accuracy - vimp::measure_accuracy(f_01, y)$point_est
acc_vim_2 <- full_accuracy - vimp::measure_accuracy(f_02, y)$point_est
acc_vim_3 <- full_accuracy - vimp::measure_accuracy(f_03, y)$point_est
acc_vim_4 <- full_accuracy - vimp::measure_accuracy(f_04, y)$point_est
acc_vim_13 <- full_accuracy - vimp::measure_accuracy(f_013, y)$point_est
acc_vim_24 <- full_accuracy - vimp::measure_accuracy(f_024, y)$point_est
acc_vims <- c(acc_vim_1, acc_vim_2, acc_vim_3, acc_vim_4, acc_vim_13, acc_vim_24)

full_auc <- cvAUC::AUC(f_0, y)
auc_vim_1 <- full_auc - cvAUC::AUC(f_01, y)
auc_vim_2 <- full_auc - cvAUC::AUC(f_02, y)
auc_vim_3 <- full_auc - cvAUC::AUC(f_03, y)
auc_vim_4 <- full_auc - cvAUC::AUC(f_04, y)
auc_vim_13 <- full_auc - cvAUC::AUC(f_013, y)
auc_vim_24 <- full_auc - cvAUC::AUC(f_024, y)
auc_vims <- c(auc_vim_1, auc_vim_2, auc_vim_3, auc_vim_4, auc_vim_13, auc_vim_24)

truths_corr <- cbind.data.frame(type = c("accuracy", "auc"),
                                  rbind(acc_vims, auc_vims))
colnames(truths_corr) <- c("type", paste0("j_", c("1", "2", "3", "4", "13_group", "24_group")))
rownames(truths_corr) <- NULL
saveRDS(truths_corr,
        here("true_vals_probit_corr.rds"))

# ------------------------------------------------------------------------------
# UNCORRELATED DATA
# ------------------------------------------------------------------------------
set.seed(1234)
x <- MASS::mvrnorm(n = sample_size, mu = x_mean, Sigma = Sigma)
y <- as.numeric((x %*% beta_0 + rnorm(sample_size)) > 0)

# compute the various population prediction functions
f_0 <- pnorm(x %*% beta_0)
f_01 <- pnorm((x[, 2] * beta_0[2]))
f_02 <- pnorm((x[, 1] * beta_0[1]))
f_03 <- f_04 <- f_0
f_013 <- pnorm(beta_0[2] * x[, 2])
f_024 <- pnorm(beta_0[1] * x[, 1])

# compute the variable importance estimates
full_accuracy <- vimp::measure_accuracy(f_0, y)$point_est
acc_vim_1 <- full_accuracy - vimp::measure_accuracy(f_01, y)$point_est
acc_vim_2 <- full_accuracy - vimp::measure_accuracy(f_02, y)$point_est
acc_vim_3 <- full_accuracy - vimp::measure_accuracy(f_03, y)$point_est
acc_vim_4 <- full_accuracy - vimp::measure_accuracy(f_04, y)$point_est
acc_vim_13 <- full_accuracy - vimp::measure_accuracy(f_013, y)$point_est
acc_vim_24 <- full_accuracy - vimp::measure_accuracy(f_024, y)$point_est
acc_vims <- c(acc_vim_1, acc_vim_2, acc_vim_3, acc_vim_4, acc_vim_13, acc_vim_24)

full_auc <- cvAUC::AUC(f_0, y)
auc_vim_1 <- full_auc - cvAUC::AUC(f_01, y)
auc_vim_2 <- full_auc - cvAUC::AUC(f_02, y)
auc_vim_3 <- full_auc - cvAUC::AUC(f_03, y)
auc_vim_4 <- full_auc - cvAUC::AUC(f_04, y)
auc_vim_13 <- full_auc - cvAUC::AUC(f_013, y)
auc_vim_24 <- full_auc - cvAUC::AUC(f_024, y)
auc_vims <- c(auc_vim_1, auc_vim_2, auc_vim_3, auc_vim_4, auc_vim_13, auc_vim_24)

truths_uncorr <- cbind.data.frame(type = c("accuracy", "auc"),
                                  rbind(acc_vims, auc_vims))
colnames(truths_uncorr) <- c("type", paste0("j_", c("1", "2", "3", "4", "13_group", "24_group")))
rownames(truths_uncorr) <- NULL

saveRDS(truths_uncorr,
        here("true_vals_probit.rds"))
