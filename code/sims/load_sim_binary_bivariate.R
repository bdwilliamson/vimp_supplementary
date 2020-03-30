#!/usr/local/bin/Rscript
## ------------------------------------------------
## FILE: load_sim_binary_bivariate.R
## CREATED: 13 March 2019 by Brian Williamson
## PURPOSE: load in simulation results
## ------------------------------------------------

## ---------------------------------------------------
## load required libraries, set up arguments
## ---------------------------------------------------
library("argparse")
library("dplyr")
library("tidyr")
library("cowplot")

parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "bivariate_loss",
                    help = "the name of the simulation")
parser$add_argument("--risk-type", default = "expected_loss",
                    help = "whether or not the risk is an expected loss")
parser$add_argument("--vimp-measure", nargs = "+", default = c("deviance", "accuracy", "auc"),
                    help = "the measure of variable importance")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 50,
                    help = "number of replicates for each job")
parser$add_argument("--b", type = "double", default = 1000,
                    help = "number of bootstrap replicates")
parser$add_argument("--cv", type = "double", default = 1,
                    help = "whether or not CV vimp estimator")
args <- parser$parse_args()

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) { # if running locally
  output_dir <- "lowdim/output/"
  plots_dir <- "lowdim/Plots/"
  code_dir <- ""
} else {
  output_dir <- "../output/"
  plots_dir <- "../Plots/"
  code_dir <- "../../"
}

## function to read with NAs
read_func <- function(x) tryCatch(readRDS(x), error = function(e) NA)
## ---------------------------------------------------
## read in the data
## ---------------------------------------------------
## list of output
output_lst <- lapply(paste0(output_dir, output_nms), read_func)
## make it a matrix
output_tib <- do.call(rbind.data.frame, output_lst)
## read in the truths
if (grepl("null", args$sim_name)) {
  truths <- readRDS(paste0(code_dir, "truths_binary_bivariate_null.rds"))
} else {
  truths <- readRDS(paste0(code_dir, "truths_binary_bivariate.rds")  )
}

if (length(args$vimp_measure) > 1) {
  if (!args$plot_all) {
    vimp_measures <- c("accuracy", "auc")
  } else {
    vimp_measures <- c("deviance")
  }
} else {
  vimp_measures <- args$vimp_measure
}
if (args$sim_name == "bivariate_naive") {
  vimp_measures <- c("deviance", "naive")
}


## ---------------------------------------------------
## compute performance metrics:
## bias*sqrt(n), var*n, mse*n, coverage
## ---------------------------------------------------
raw_performance <- output_tib %>%
  filter(type %in% vimp_measures) %>%
  mutate(bias = (est - truth)*sqrt(n),
         mse = (est - truth)^2*n, cover = (cil <= truth & ciu >= truth))

## average over everything
average_performance <- raw_performance %>%
  group_by(n, j, type, delta) %>%
  # select(-truth, -est, -se, -cil, -ciu, -mc_id) %>%
  select(-truth, -se, -cil, -ciu, -risk_full, -risk_reduced) %>%
  summarize(bias = mean(bias, na.rm = TRUE), var = var(est, na.rm = TRUE),
            mse = mean(mse, na.rm = TRUE), cover = mean(cover, na.rm = TRUE),
            power = mean(test, na.rm = TRUE),
            sd = sd(est, na.rm = TRUE)) %>%
  ungroup()
if (args$sim_name == "bivariate_naive") {
  average_performance <- average_performance %>%
    group_by(n, j, type, delta) %>%
    filter(!is.na(n)) %>%
    mutate(sd = sd*sqrt(n)/sqrt(args$nreps_total), var = var*n) %>%
    ungroup()
} else {
  average_performance <- average_performance %>%
    filter(n != 300) %>%
    group_by(n, j, type, delta) %>%
    filter(!is.na(n)) %>%
    mutate(sd = sd*sqrt(n)/sqrt(args$nreps_total), var = var*n) %>%
    ungroup()

}

## ---------------------------------------------------
## make plots
## ---------------------------------------------------
point_size <- 5
text_size <- 20
title_text_size <- 30
y_lim_bias <- c(-0.35, 0.35)
y_lim_mse <- c(0, 1)
y_lim_var <- c(0, 2)
if (grepl("naive", args$sim_name)) {
  dodge_x <- 100
} else {
  dodge_x <- 300
}
right_pad <- 10
y_lim_cover <- c(0.5, 1)
legend_pos <- c(0.1, 0.8)
if (args$sim_name == "bivariate_naive") {
  point_vals <- c(16, 16, 15, 15)
  color_vals <- rep(c("black", "blue"), 2)
} else {
  point_vals <- rep(c(16, 15, 1, 4, 18, 8))[1:(2*length(vimp_measures))]
  # point_vals <-  c(16, 16, 15, 15, 18, 18)[1:(2*length(vimp_measures))]
  color_vals <- rep(c("black", "blue", "red")[1:ifelse(length(vimp_measures) == 1, 2, length(vimp_measures))], each = length(vimp_measures))
}
if (grepl("null", args$sim_name) | length(vimp_measures) == 1) {
  y_lim_mse <- c(0, 2.2)
}
if (args$sim_name == "bivariate_naive" | args$cv == 0 | grepl("simple_lib", args$sim_name)) {
  y_lim_mse <- c(0, 2)
}
if (grepl("null", args$sim_name) & !grepl("simple_lib", args$sim_name) & args$cv == 0) {
  y_lim_cover <- c(0, 1)
  y_lim_mse <- c(0, 200)
}
if (args$plot_all) {
  y_lim_mse <- c(0, 5)
}

bias_plot <- average_performance %>%
  filter(delta == 0) %>%
  ggplot(aes(x = n, y = bias, group = factor(paste(n, j, type, sep = "_")),
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  geom_errorbar(aes(ymin = bias - 1.96*sd, ymax = bias + 1.96*sd), width = rep(200, dim(average_performance)[1]),
                position = position_dodge(width = dodge_x), size = 0.3*point_size) +
  ylim(y_lim_bias) +
  ylab(expression(paste(sqrt(n), "x empirical ", bias[n]))) +
  xlab("n") +
  ggtitle(expression(bold(paste("Empirical bias scaled by ", sqrt(n))))) +
  scale_color_manual(name = "Measure; feature", values = color_vals) +
  scale_shape_manual(name = "Measure; feature", values = point_vals) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  guides(color = FALSE, shape = FALSE) +
  theme(text = element_text(size = text_size),
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

variance_plot <- average_performance %>%
  filter(delta == 0) %>%
  ggplot(aes(x = n, y = var, group = factor(paste(n, j, type, sep = "_")),
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab(expression(paste(n, "x empirical ", var[n]))) +
  ylim(y_lim_var) +
  xlab("n") +
  ggtitle(expression(bold(paste("Empirical variance scaled by ", n)))) +
  scale_color_manual(name = "Measure; feature", values = color_vals) +
  scale_shape_manual(name = "Measure; feature", values = point_vals) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  guides(color = FALSE, shape = FALSE) +
  theme(text = element_text(size = text_size),
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

mse_plot <- average_performance %>%
  filter(delta == 0) %>%
  ggplot(aes(x = n, y = mse, group = factor(paste(n, j, type, sep = "_")),
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab(expression(paste(n, "x empirical ", MSE[n]))) +
  ylim(y_lim_mse) +
  xlab("n") +
  scale_color_manual(name = "Measure; feature", values = color_vals) +
  scale_shape_manual(name = "Measure; feature", values = point_vals) +
  ggtitle(expression(bold(paste("Empirical mean squared error scaled by ", n)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(legend.position = legend_pos, legend.box.background = element_rect(colour = "black"),
        text = element_text(size = text_size),
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))


cover_plot <- average_performance %>%
  filter(delta == 0) %>%
  ggplot(aes(x = n, y = cover, group = factor(paste(n, j, type, sep = "_")),
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab("Empirical coverage") +
  xlab("n") +
  scale_color_manual(name = "Measure; feature", values = color_vals) +
    scale_shape_manual(name = "Measure; feature", values = point_vals) +
    ggtitle("Empirical coverage of nominal 95% CIs") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(y_lim_cover) +
  guides(color = FALSE, shape = FALSE) +
  theme(text = element_text(size = text_size),
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

power_plot <- average_performance %>%
  filter(delta == 0) %>%
  ggplot(aes(x = n, y = power, group = factor(paste(n, j, type, sep = "_")),
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab("Empirical power") +
  xlab("n") +
  scale_color_manual(name = "Measure; feature", values = color_vals) +
  scale_shape_manual(name = "Measure; feature", values = point_vals) +
  ggtitle(expression(bold(paste("Proportion of tests rejected (", beta, " = ", 0, ")", sep = "")))) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  ylim(c(0, 1)) +
  guides(color = FALSE, shape = FALSE) +
  theme(text = element_text(size = text_size),
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

if (length(args$delta) > 1) {
  power_plot_2 <- average_performance %>%
    filter(delta == args$delta[2]) %>%
    ggplot(aes(x = n, y = power, group = factor(paste(n, j, type, sep = "_")),
               shape = factor(paste(type, j, sep = "; ")),
               color = factor(paste(type, j, sep = "; ")))) +
    geom_point(position = position_dodge(width = dodge_x), size = point_size) +
    ylab("Empirical power") +
    xlab("n") +
    scale_color_manual(name = "Measure; feature", values = color_vals) +
    scale_shape_manual(name = "Measure; feature", values = point_vals) +
    ggtitle(expression(bold(paste("Proportion of tests rejected (", beta, " = ", 0.05, ")", sep = "")))) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    ylim(c(0, 1)) +
    guides(color = FALSE, shape = FALSE) +
    theme(text = element_text(size = text_size),
          axis.text = element_text(size = text_size),
          plot.title = element_text(size = text_size),
          plot.margin = unit(c(0, right_pad, 0, 0), "mm"))
}

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  if (length(args$delta) == 1) {
    plot_grid(mse_plot, cover_plot, power_plot)
  } else {
    plot_grid(mse_plot, cover_plot, power_plot, power_plot_2)
  }
  if (grepl("naive", args$sim_name)) {
    plot_grid(bias_plot, cover_plot, power_plot)
  }
}

fig_width <- 15
fig_height <- 0.75*fig_width
if (grepl("naive", args$sim_name)) {
  ggsave(filename = paste0(plots_dir, args$sim_name, "_performance_", paste0(vimp_measures, collapse = "_"), "_cv_", args$cv, ".png"),
         plot = plot_grid(bias_plot, cover_plot, variance_plot, power_plot),
         width = fig_width, height = fig_height)
  ggsave(filename = paste0(plots_dir, args$sim_name, "_performance_", paste0(vimp_measures, collapse = "_"), "_cv_", args$cv, ".tiff"),
         plot = plot_grid(bias_plot, cover_plot, variance_plot, power_plot),
         width = fig_width, height = fig_height)
} else {
  if (length(args$delta) > 1) {
    plt <- plot_grid(mse_plot, cover_plot, power_plot, power_plot_2)
  } else {
    plt <- plot_grid(mse_plot, cover_plot, power_plot)
  }
  ggsave(filename = paste0(plots_dir, args$sim_name, "_performance_", paste0(vimp_measures, collapse = "_"), "_cv_", args$cv, ".png"),
         plot = plt,
         width = fig_width, height = fig_height)
  ggsave(filename = paste0(plots_dir, args$sim_name, "_performance_", paste0(vimp_measures, collapse = "_"), "_cv_", args$cv, ".tiff"),
         plot = plt,
         width = fig_width, height = fig_height)

}

## -----------------------------------------------------
## assess trends in estimates
## -----------------------------------------------------

overall_means <- output_tib %>%
  group_by(n, j, type) %>%
  summarize(truth = mean(truth), est = mean(est), se = mean(se),
            cil = mean(cil), ciu = mean(ciu), risk_full = mean(risk_full), risk_reduced = mean(risk_reduced))
overall_means
