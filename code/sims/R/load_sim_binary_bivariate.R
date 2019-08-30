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
## names of files to read in
if (args$cv == 1) {
  output_nms <- paste0(args$sim_name, "_", paste0(args$vimp_measure, collapse = "_"), "_", 1:(args$nreps_total/args$nreps_per_job * 20), ".rds")  
} else {
  output_nms <- paste0("_no_cv_results/", args$sim_name, "_", paste0(args$vimp_measure, collapse = "_"), "_cv_", args$cv, "_", 1:(args$nreps_total/args$nreps_per_job * 20), ".rds")
}

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
## make it tidy if I have naive
if (args$sim_name == "bivariate_naive") {
  output_naive <- output_tib %>% 
    select(n, j, truth, mc_id, naive_deviance, naive_cil_deviance, naive_ciu_deviance) %>% 
    mutate(type = "naive", est = naive_deviance, cil = naive_cil_deviance, ciu = naive_ciu_deviance,
           se = NA, test = NA, p_value = NA, risk_full = NA, risk_reduced = NA) %>% 
    select(-naive_deviance, -naive_cil_deviance, -naive_ciu_deviance)
  output_risk <- output_tib %>% 
    select(-naive_deviance, -naive_cil_deviance, -naive_ciu_deviance) %>% 
    mutate(type = "predictiveness")
  output_tib_2 <- bind_rows(output_risk, output_naive)
  output_tib <- output_tib_2 
}

  
## ---------------------------------------------------
## compute performance metrics:
## bias*sqrt(n), var*n, mse*n, coverage
## ---------------------------------------------------
raw_performance <- output_tib %>% 
  mutate(bias = (est - truth)*sqrt(n),
         mse = (est - truth)^2*n, cover = (cil <= truth & ciu >= truth))
  
## average over everything 
average_performance <- raw_performance %>% 
  group_by(n, j, type) %>% 
  # select(-truth, -est, -se, -cil, -ciu, -mc_id) %>% 
  select(-truth, -se, -cil, -ciu, -risk_full, -risk_reduced) %>% 
  summarize(bias = mean(bias, na.rm = TRUE), var = var(est, na.rm = TRUE),
            mse = mean(mse, na.rm = TRUE), cover = mean(cover, na.rm = TRUE),
            power = mean(test, na.rm = TRUE),
            sd = sd(est, na.rm = TRUE)) %>% 
  ungroup()
average_performance <- average_performance %>% 
  filter(n != 300) %>% 
  group_by(n, j, type) %>% 
  filter(!is.na(n)) %>% 
  mutate(sd = sd*sqrt(n)/sqrt(args$nreps_total), var = var*n) %>% 
  ungroup()

## ---------------------------------------------------
## make plots
## ---------------------------------------------------
point_size <- 5
text_size <- 20
title_text_size <- 30
y_lim_bias <- c(-0.35, 0.35)
y_lim_mse <- c(0, 2)
y_lim_var <- c(0, 2)
dodge_x <- 300
dodge_x_large <- 400
right_pad <- 10
if (args$sim_name == "bivariate_naive") {
  point_vals <- c(16, 13, 18, 9)
  color_vals <- rep(c("black", "blue"), 2)
  legend_pos <- c(0.1, 0.325)
} else {
  point_vals <-  c(16, 13, 18, 9, 15, 7)[1:(2*length(args$vimp_measure))]
  color_vals <- rep(c("black", "blue"), length(args$vimp_measure))
  legend_pos <- c(0.7, 0.325)
}

bias_plot <- average_performance %>% 
  ggplot(aes(x = n, y = bias, group = factor(paste(n, j, type, sep = "_")), 
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  geom_errorbar(aes(ymin = bias - 1.96*sd, ymax = bias + 1.96*sd), width = rep(200, dim(average_performance)[1]),
                position = position_dodge(width = dodge_x), size = 0.3*point_size) +
  ylim(y_lim_bias) +
  ylab(expression(paste(sqrt(n), "x estimated ", bias[n]))) +
  xlab("n") +
  ggtitle(expression(bold(paste("Estimated bias scaled by ", sqrt(n))))) +
  scale_color_manual(name = "Measure; j", values = color_vals) +
  scale_shape_manual(name = "Measure; j", values = point_vals) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  guides(color = FALSE, shape = FALSE) +
  theme(text = element_text(size = text_size), 
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

variance_plot <- average_performance %>% 
  ggplot(aes(x = n, y = var, group = factor(paste(n, j, type, sep = "_")), 
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab(expression(paste(n, "x estimated ", var[n]))) +
  ylim(y_lim_var) +
  xlab("n") +
  ggtitle(expression(bold(paste("Estimated variance scaled by ", n)))) +
  scale_color_manual(name = "Measure; j", values = color_vals) +
  scale_shape_manual(name = "Measure; j", values = point_vals) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  guides(color = FALSE, shape = FALSE) +
  theme(text = element_text(size = text_size), 
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))
  
mse_plot <- average_performance %>% 
  ggplot(aes(x = n, y = mse, group = factor(paste(n, j, type, sep = "_")), 
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab(expression(paste(n, "x estimated ", MSE[n]))) +
  ylim(y_lim_mse) +
  xlab("n") +
  scale_color_manual(name = "Measure; j", values = color_vals) +
  scale_shape_manual(name = "Measure; j", values = point_vals) +
  ggtitle(expression(bold(paste("Estimated mean squared error scaled by ", n)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme(text = element_text(size = text_size), 
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

cover_plot <- average_performance %>% 
  ggplot(aes(x = n, y = cover, group = factor(paste(n, j, type, sep = "_")), 
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab("Coverage") +
  xlab("n") +
  scale_color_manual(name = "Measure; j", values = color_vals) +
    scale_shape_manual(name = "Measure; j", values = point_vals) +
    ggtitle("Coverage of nominal 95% CIs") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  ylim(c(0, 1)) +
  theme(legend.position = legend_pos, legend.box.background = element_rect(colour = "black"),
        text = element_text(size = text_size), 
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        legend.text = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

power_plot <- average_performance %>% 
  ggplot(aes(x = n, y = power, group = factor(paste(n, j, type, sep = "_")), 
             shape = factor(paste(type, j, sep = "; ")),
             color = factor(paste(type, j, sep = "; ")))) +
  geom_point(position = position_dodge(width = dodge_x), size = point_size) +
  ylab("Power") +
  xlab("n") +
  scale_color_manual(name = "Measure; j", values = color_vals) +
  scale_shape_manual(name = "Measure; j", values = point_vals) +
  ggtitle("Proportion of tests rejected") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  ylim(c(0, 1)) +
  guides(color = FALSE, shape = FALSE) +
  theme(text = element_text(size = text_size), 
        axis.text = element_text(size = text_size),
        plot.title = element_text(size = text_size),
        plot.margin = unit(c(0, right_pad, 0, 0), "mm"))

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  plot_grid(bias_plot, cover_plot, variance_plot, power_plot)
  # mse_plot
}

fig_width <- 15
fig_height <- 0.75*fig_width
ggsave(filename = paste0(plots_dir, args$sim_name, "_performance_", paste0(args$vimp_measure, collapse = "_"), "_cv_", args$cv, ".png"),
       plot = plot_grid(bias_plot, cover_plot, variance_plot, power_plot),
       width = fig_width, height = fig_height)
ggsave(filename = paste0(plots_dir, args$sim_name, "_performance_", paste0(args$vimp_measure, collapse = "_"), "_cv_", args$cv, ".tiff"),
       plot = plot_grid(bias_plot, cover_plot, variance_plot, power_plot),
       width = fig_width, height = fig_height)


## -----------------------------------------------------
## assess trends in estimates
## -----------------------------------------------------

overall_means <- output_tib %>% 
  group_by(n, j, type) %>% 
  summarize(truth = mean(truth), est = mean(est), se = mean(se), 
            cil = mean(cil), ciu = mean(ciu), risk_full = mean(risk_full), risk_reduced = mean(risk_reduced))
overall_means  
