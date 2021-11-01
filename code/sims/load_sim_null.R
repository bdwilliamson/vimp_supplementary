# load the simulation investigating null properties in low dimensions

#### load required libraries and set up arguments ####
library("argparse")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("tibble")
library("here")
library("data.table")

code_dir <- "sim_code"
source(here(code_dir, "utils.R"))

parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "investigate-null-props",
                    help = "the name of the simulation")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 5,
                    help = "number of replicates for each job")
parser$add_argument("--bootstrap", default = 0, type = "double",
                    help = "should we use the bootstrap for variance estimation?")
parser$add_argument("--p", default = 4, type = "double",
                    help = "dimension")
parser$add_argument("--use-cf-se", default = 1, type = "double",
                    help = "should we use the cross-fitted SE estimator?")
args <- parser$parse_args()

all_estimators <- c("glm", "gam", "ranger", "SL")
all_cvs <- c(0, 1)
all_js <- c(3, 2)
all_vims <- c("accuracy", "auc")
main_vim <- "accuracy"
main_j <- 3
#### read in the results ####
nms_grid <- expand.grid(est = all_estimators, cv = all_cvs, j = all_js)
output_nms <- here("sim_output", args$sim_name, 
                   paste0(nms_grid$est, "_cv_", nms_grid$cv, "_j_", nms_grid$j,
                          "_p_4_corr_0_boot_", args$bootstrap, ".rds"))
output_tib <- tibble::tibble(data.table::rbindlist(
  lapply(as.list(output_nms), readRDS), fill = TRUE
))
args$corr <- grepl("corr", args$sim_name)
truth <- get_true_values(args, num_digits = 3, dir = here(code_dir)) %>% 
  pivot_longer(cols = starts_with("j_"), names_to = "j",
               names_prefix = "j_")
if (args$p == 4) {
  truth <- truth %>% 
    filter(j %in% c("1", "2", "3"))
}
output_tib_fixed_truths <- output_tib %>% 
  select(-truth) %>% 
  left_join(truth %>% rename(truth = value), by = c("type", "j"))
# if we aren't using the cross-fitted SE estimator, need to recompute CIs and p-values
if (!as.logical(args$use_cf_se)) {
  output_tib_fixed_truths <- output_tib_fixed_truths %>% 
    mutate(cil = ifelse(cv == 1, est - qnorm(0.975) * non_cf_se, cil),
           ciu = ifelse(cv == 1, est + qnorm(0.975) * non_cf_se, ciu),
           test = ifelse(cv == 1, 1 - pnorm(est / non_cf_se) < 0.05, test))
}
#### summarize performance ####
performance_tib <- output_tib_fixed_truths %>% 
  mutate(bias_init = sqrt(n) * (est - truth),
         mse_init = n * (est - truth) ^ 2,
         cover_init = cil <= truth & ciu >= truth,
         type_1_error_init = ifelse(se == 0, FALSE, test),
         width_init = ciu - cil) %>% 
  group_by(n, j, estimator, cv, type) %>% 
  summarize(bias = mean(bias_init), mse = mean(mse_init), 
            cover = mean(cover_init), width = mean(width_init),
            type_1_error = mean(type_1_error_init, na.rm = TRUE),
            var_init = var(est), 
            bias_mc_se = sqrt( mean((bias_init - bias) ^ 2) / args$nreps_total), 
            .groups = "drop") %>% 
  mutate(variance = n * var_init) %>% 
  select(-var_init)

plot_tib <- performance_tib %>% 
  mutate(estimator = factor(estimator, labels = c("GAM", "GLM", "RF", "SL"),
                            levels = c("gam", "glm", "ranger", "SL")), 
         n = factor(n), 
         cv = factor(cv, labels = c("Not cross-fitted", "Cross-fitted"), 
                     levels = c(0, 1)))

#### create individual panels for the final plots ####
unique_ns <- unique(performance_tib$n)
point_size <- 1.5
title_text_size <- 12
legend_text_size <- 10
axis_text_size <- 10
fig_width <- 9
fig_height <- 4.25
dodge_width <- 0.875
color_begin <- 0
color_end <- 0.5
label_x <- 0.025
all_mins_and_maxes <- plot_tib %>%
  summarize(across(c(bias, width, variance),
                   .fns = list(max = max, min = min)))
bias_lim <- round(max(abs(all_mins_and_maxes$bias_min), abs(all_mins_and_maxes$bias_max)), 2)
width_lim <- c(0, all_mins_and_maxes$width_max)
variance_lim <- c(0, all_mins_and_maxes$variance_max)
for (this_j in all_js) {
  test_y_min <- ifelse(this_j == 3, 0, 
                       round(min(plot_tib %>% 
                                   filter(j == this_j) %>% 
                                   pull(type_1_error)) - 0.1, 2))
  test_y_max <- ifelse(this_j != 3, 1, 
                       round(max(plot_tib %>% 
                                   filter(j == this_j) %>% 
                                   pull(type_1_error)) + 0.1, 2))
  for (vim in all_vims) {
    in_main_ms <- (vim == main_vim) & (this_j == main_j) & (args$bootstrap == 0) &
      !(args$use_cf_se == 0)
    if (!in_main_ms) {
      # make common y-axis limits for each plot
      this_bias_lim <- bias_lim
      this_var_lim <- variance_lim
    } else {
      this_bias_lim <- plot_tib %>% filter(type == vim, j == this_j) %>%
        summarize(bias = max(abs(bias))) %>% pull(bias)
      this_var_lim_upper <- plot_tib %>% filter(type == vim, j == this_j) %>%
        summarize(var = max(abs(variance))) %>% pull(var)
      this_var_lim <- c(0, this_var_lim_upper)
    }
    bias_plot_j <- plot_tib %>%
      filter(type == vim, j == this_j) %>%
      ggplot(aes(x = estimator, y = bias, shape = estimator, color = cv)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_point(position = position_dodge(dodge_width), size = point_size) +
      scale_color_viridis_d(begin = 0, end = 0.5) +
      ylim(c(-this_bias_lim, this_bias_lim)) +
      ylab(expression(paste(sqrt(n), " x empirical ", bias[n]))) +
      xlab("n") +
      labs(shape = "Estimator", color = "Cross-fitting") +
      ggtitle("BIAS") +
      facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
            panel.spacing = unit(0, "cm"),
            strip.background = element_blank(), strip.placement = "outside",
            panel.grid.minor.x = element_line(color = "grey85"),
            panel.grid.major.y = element_line(color = "grey85")) +
      geom_vline(aes(xintercept = 0.4), color = "grey85")
    var_plot_j <- plot_tib %>%
      filter(type == vim, j == this_j) %>%
      ggplot(aes(x = estimator, y = variance, shape = estimator, color = cv)) +
      geom_point(position = position_dodge(dodge_width), size = point_size) +
      scale_color_viridis_d(begin = 0, end = 0.5) +
      ylim(this_var_lim) +
      ylab(expression(paste(n, " x empirical ", variance[n]))) +
      xlab("n") +
      labs(shape = "Estimator", color = "Cross-fitting") +
      ggtitle("VARIANCE") +
      facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
            panel.spacing = unit(0, "cm"),
            strip.background = element_blank(), strip.placement = "outside",
            panel.grid.minor.x = element_line(color = "grey85"),
            panel.grid.major.y = element_line(color = "grey85")) +
      geom_vline(aes(xintercept = 0.4), color = "grey85")
    cover_plot_j <- plot_tib %>%
      filter(type == vim, j == this_j) %>%
      ggplot(aes(x = estimator, y = cover, shape = estimator, color = cv)) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
      geom_point(position = position_dodge(dodge_width), size = point_size) +
      scale_color_viridis_d(begin = 0, end = 0.5) +
      ylim(c(0.5, 1)) +
      ylab("Empirical coverage") +
      xlab("n") +
      labs(shape = "Estimator", color = "Cross-fitting") +
      ggtitle("COVERAGE") +
      facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
            panel.spacing = unit(0, "cm"),
            strip.background = element_blank(), strip.placement = "outside",
            panel.grid.minor.x = element_line(color = "grey85"),
            panel.grid.major.y = element_line(color = "grey85")) +
      geom_vline(aes(xintercept = 0.4), color = "grey85")
    testing_intercept <- switch((this_j == 2) + 1, 0.05, 1)
    testing_title <- ifelse(this_j == 3, "TYPE I ERROR", "POWER")
    alpha_plot_j <- plot_tib %>% 
      filter(type == vim, j == this_j) %>% 
      ggplot(aes(x = estimator, y = type_1_error, shape = estimator, color = cv)) +
      geom_hline(yintercept = testing_intercept, linetype = "dashed", color = "red") +
      geom_point(position = position_dodge(dodge_width), size = point_size) +
      scale_color_viridis_d(begin = color_begin, end = color_end) +
      ylim(c(test_y_min, test_y_max)) +
      ylab("Proportion of tests rejected") +
      xlab("n") +
      labs(shape = "Estimator", color = "Cross-fitting") +
      ggtitle(testing_title) +
      facet_wrap(~ n, nrow = 1, labeller = "label_value", strip.position = "bottom") +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
            panel.spacing = unit(0, "cm"),
            strip.background = element_blank(), strip.placement = "outside",
            panel.grid.minor.x = element_line(color = "grey85"),
            panel.grid.major.y = element_line(color = "grey85")) +
      geom_vline(aes(xintercept = 0.4), color = "grey85")
    four_panel_plot_j <- plot_grid(
      bias_plot_j + theme(legend.position = "none", 
                          title = element_text(size = title_text_size),
                          axis.title = element_text(size = axis_text_size),
                          axis.text = element_text(size = axis_text_size),
                          plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
      var_plot_j + theme(legend.position = "none", 
                         title = element_text(size = title_text_size),
                         axis.title = element_text(size = axis_text_size),
                         axis.text = element_text(size = axis_text_size),
                         plot.margin = unit(c(0, 0, 0, 0), "cm")),
      cover_plot_j + theme(legend.position = "none", title = element_text(size = title_text_size),
                           axis.title = element_text(size = axis_text_size),
                           axis.text = element_text(size = axis_text_size),
                           plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
      alpha_plot_j + theme(legend.position = "none", title = element_text(size = title_text_size),
                           axis.title = element_text(size = axis_text_size),
                           axis.text = element_text(size = axis_text_size),
                           plot.margin = unit(c(0, 0, 0, 0), "cm")),
      labels = "AUTO", label_x = label_x, label_size = title_text_size
    )
    if (args$use_cf_se == 0) {
      four_panel_plot_j <- plot_grid(cover_plot_j + theme(legend.position = "none", 
                                                          title = element_text(size = title_text_size),
                                                          axis.title = element_text(size = axis_text_size),
                                                          axis.text = element_text(size = axis_text_size),
                                                          plot.margin = unit(c(0, 0, 0, 0), "cm")),
                                     alpha_plot_j + theme(legend.position = "none", 
                                                          title = element_text(size = title_text_size),
                                                          axis.title = element_text(size = axis_text_size),
                                                          axis.text = element_text(size = axis_text_size),
                                                          plot.margin = unit(c(0, 0, 0, 0), "cm")),
                                     labels = "AUTO", label_x = label_x)
    }
    legend_j <- get_legend(
      bias_plot_j + 
        guides(shape = guide_legend(nrow = 2),
               color = guide_legend(nrow = 2)) +
        theme(legend.direction = "horizontal",
              legend.position = "bottom",
              legend.title = element_text(size = legend_text_size),
              legend.text = element_text(size = legend_text_size),
              legend.spacing.x = unit(1.5, "cm"))
    )
    full_plot_j <- plot_grid(four_panel_plot_j, legend_j, ncol = 1, nrow = 2,
                             rel_heights = c(1, .1))
    plot_name <- paste0(
      switch((vim == main_vim & this_j == main_j & !(args$use_cf_se == 0)) + 1, 
             "supplement_", "main_"),
      "null_", vim, "_", this_j, "_p_4_boot_", args$bootstrap, 
      switch((args$use_cf_se == 0) + 1, "", "_nocf")
    )
    ggsave(filename = here("plots", paste0(plot_name, ".png")),
           plot = full_plot_j, device = "png", width = fig_width, height = fig_height, dpi = 300)
    ggsave(filename = here("plots", paste0(plot_name, ".eps")),
           plot = full_plot_j, device = "eps", width = fig_width, height = fig_height, dpi = 300)
    if (!as.logical(args$use_cf_se)) {
      saveRDS(four_panel_plot_j, here("plots", paste0(plot_name, ".rds")))
    }
  }
}
