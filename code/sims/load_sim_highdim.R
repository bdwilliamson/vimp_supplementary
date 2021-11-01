#!/usr/local/bin/Rscript
# -------------------------------------------------------
# load the results from the higher-dimensional simulation
# -------------------------------------------------------

#### load required libraries and set up arguments ####
library("argparse")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("data.table")
library("tibble")
library("here")

parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "investigate-highdim-props",
                    help = "the name of the simulation")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates in total")
parser$add_argument("--nreps-per-job", type = "double", default = 2,
                    help = "number of replicates for each job")
args <- parser$parse_args()

all_estimators <- c("SL")
all_cvs <- c(1)
all_corrs <- c(0, 1)
all_features <- c("1", "2", "3", "4", "1,3", "2,4")
all_ps <- c(50, 100, 200)
all_vims <- c("accuracy", "auc")

#### read in the results ####
nms_grid <- expand.grid(est = all_estimators, cv = all_cvs, 
                        p = all_ps, corr = all_corrs, s = all_features)
output_nms <- here("sim_output", args$sim_name, 
                   paste0(nms_grid$est, "_cv_", nms_grid$cv, "_j_", nms_grid$s,
                          "_p_", nms_grid$p, "_corr_", nms_grid$corr, "_boot_0.rds"))
output_tib <- tibble::tibble(data.table::rbindlist(
  lapply(as.list(output_nms), readRDS), fill = TRUE
))

#### summarize performance ####
performance_tib <- output_tib %>% 
  mutate(bias_init = sqrt(n) * (est - truth),
         mse_init = n * (est - truth) ^ 2,
         cover_init = cil <= truth & ciu >= truth,
         type_1_error_init = ifelse(se == 0, FALSE, test),
         width_init = ciu - cil) %>% 
  group_by(n, p, corr, j, estimator, cv, type) %>% 
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
         p = factor(p, levels = c(50, 100, 200), ordered = TRUE),
         corr_fct = factor(corr, labels = c("No correlation", "Correlation"),
                       levels = c(0, 1)),
         legend_ord = case_when(
           j == "1" ~ 1,
           j == "2" ~ 2,
           j == "3" ~ 5,
           j == "4" ~ 6,
           j == "1,3" ~ 3,
           j == "2,4" ~ 4
         ), 
         group_fct = factor(paste(p, j, sep = "_")),
         shape_fct = factor(j))

#### create individual panels for the final plots ####
unique_ps <- unique(performance_tib$p)
point_size <- 1.5
title_text_size <- 12
legend_text_size <- 12
axis_text_size <- 12
fig_width <- 9
fig_height <- 7
dodge_width <- 0.875
color_begin <- 0
color_end <- 0.5
label_x <- 0.025
all_mins_and_maxes <- plot_tib %>%
  summarize(across(c(bias, width, variance, cover),
                   .fns = list(max = max, min = min)))
bias_lim <- round(max(abs(all_mins_and_maxes$bias_min), abs(all_mins_and_maxes$bias_max)), 2)
width_lim <- c(0, all_mins_and_maxes$width_max)
variance_lim <- c(0, all_mins_and_maxes$variance_max)
cover_lim <- round(c(abs(all_mins_and_maxes$cover_min), 
                       max(abs(all_mins_and_maxes$cover_max), 1)), 2)
for (this_corr in all_corrs) {
  for (vim in all_vims) {
    this_bias_plot <- plot_tib %>% 
      filter(type == vim, corr == this_corr) %>% 
      ggplot(aes(x = shape_fct, y = bias, shape = shape_fct, group = group_fct)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_point(position = position_dodge(dodge_width), size = point_size) +
      scale_color_viridis_d(begin = color_begin, end = color_end) +
      ylim(c(-bias_lim, bias_lim)) +
      ylab(expression(paste(sqrt(n), " x empirical ", bias[n]))) +
      labs(x = NULL) +
      labs(shape = "Feature(s)") +
      ggtitle("BIAS") +
      facet_grid(rows = vars(p), cols = vars(n), 
                 labeller = label_both, scales = "free",
                 switch = "both") +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0.2, "cm"),
            strip.background = element_blank(), strip.placement = "outside",
            panel.grid.minor.x = element_line(color = "grey85"),
            panel.grid.major.y = element_line(color = "grey85")) +
      geom_vline(aes(xintercept = 0.4), color = "grey85")
    this_var_plot <- plot_tib %>% 
      filter(type == vim, corr == this_corr) %>% 
      ggplot(aes(x = shape_fct, y = variance, shape = shape_fct, group = group_fct)) +
      geom_point(position = position_dodge(dodge_width), size = point_size) +
      scale_color_viridis_d(begin = color_begin, end = color_end) +
      ylim(variance_lim) +
      ylab(expression(paste(n, " x empirical ", variance[n]))) +
      labs(x = NULL) +
      labs(shape = "Feature(s)") +
      ggtitle("VARIANCE") +
      facet_grid(rows = vars(p), cols = vars(n), 
                 labeller = label_both, scales = "free",
                 switch = "both") +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0.2, "cm"),
            strip.background = element_blank(), strip.placement = "outside",
            panel.grid.minor.x = element_line(color = "grey85"),
            panel.grid.major.y = element_line(color = "grey85")) +
      geom_vline(aes(xintercept = 0.4), color = "grey85")
    this_cover_plot <- plot_tib %>% 
      filter(type == vim, corr == this_corr) %>% 
      ggplot(aes(x = shape_fct, y = cover, shape = shape_fct, group = group_fct)) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
      geom_point(position = position_dodge(dodge_width), size = point_size) +
      scale_color_viridis_d(begin = color_begin, end = color_end) +
      ylim(cover_lim) +
      ylab("Empirical coverage") +
      labs(x = NULL) +
      labs(shape = "Feature(s)") +
      ggtitle("COVERAGE") +
      facet_grid(rows = vars(p), cols = vars(n), 
                 labeller = label_both, scales = "free",
                 switch = "both") +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0.2, "cm"),
            strip.background = element_blank(), strip.placement = "outside",
            panel.grid.minor.x = element_line(color = "grey85"),
            panel.grid.major.y = element_line(color = "grey85")) +
      geom_vline(aes(xintercept = 0.4), color = "grey85")
    if (!grepl("alt", args$sim_name)) {
      this_fourth_plot <- plot_tib %>% 
        filter(type == vim, corr == this_corr) %>% 
        ggplot(aes(x = shape_fct, y = type_1_error, shape = shape_fct, group = group_fct)) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
        geom_point(position = position_dodge(dodge_width), size = point_size) +
        scale_color_viridis_d(begin = color_begin, end = color_end) +
        ylab("Proportion of tests rejected") +
        ylim(c(0, 1)) +
        labs(x = NULL) +
        labs(shape = "Feature(s)") +
        ggtitle("TESTING") +
        facet_grid(rows = vars(p), cols = vars(n), 
                   labeller = label_both, scales = "free",
                   switch = "both") +
        theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0.2, "cm"),
              strip.background = element_blank(), strip.placement = "outside",
              panel.grid.minor.x = element_line(color = "grey85"),
              panel.grid.major.y = element_line(color = "grey85")) +
        geom_vline(aes(xintercept = 0.4), color = "grey85")  
    } else {
      this_fourth_plot <- plot_tib %>% 
        filter(type == vim, corr == this_corr) %>% 
        ggplot(aes(x = shape_fct, y = width, shape = shape_fct, group = group_fct)) +
        geom_point(position = position_dodge(dodge_width), size = point_size) +
        scale_color_viridis_d(begin = color_begin, end = color_end) +
        ylim(width_lim) +
        ylab("Confidence interval width") +
        labs(x = NULL) +
        labs(shape = "Feature(s)") +
        ggtitle("WIDTH") +
        facet_grid(rows = vars(p), cols = vars(n), 
                   labeller = label_both, scales = "free",
                   switch = "both") +
        theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0, "cm"),
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0.2, "cm"),
              strip.background = element_blank(), strip.placement = "outside",
              panel.grid.minor.x = element_line(color = "grey85"),
              panel.grid.major.y = element_line(color = "grey85")) +
        geom_vline(aes(xintercept = 0.4), color = "grey85")  
    }
    
    this_four_panel_plot <- plot_grid(
      this_bias_plot + theme(legend.position = "none", 
                             title = element_text(size = title_text_size),
                             axis.title = element_text(size = axis_text_size),
                             axis.text = element_text(size = axis_text_size),
                             plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
      this_var_plot + theme(legend.position = "none", 
                            title = element_text(size = title_text_size),
                            axis.title = element_text(size = axis_text_size),
                            axis.text = element_text(size = axis_text_size),
                            plot.margin = unit(c(0, 0, 0, 0), "cm")),
      this_cover_plot + theme(legend.position = "none", 
                              title = element_text(size = title_text_size),
                              axis.title = element_text(size = axis_text_size),
                              axis.text = element_text(size = axis_text_size),
                              plot.margin = unit(c(0.1, 0, 0, 0), "cm")),
      this_fourth_plot + theme(legend.position = "none", 
                               title = element_text(size = title_text_size),
                               axis.title = element_text(size = axis_text_size),
                               axis.text = element_text(size = axis_text_size),
                               plot.margin = unit(c(0, 0, 0, 0), "cm")),
      labels = "AUTO", label_x = label_x, label_size = title_text_size
    )
    this_legend <- get_legend(
      this_bias_plot + 
        guides(shape = guide_legend(nrow = 2),
               color = guide_legend(nrow = 2)) +
        theme(legend.direction = "horizontal",
              legend.position = "bottom",
              legend.title = element_text(size = legend_text_size),
              legend.text = element_text(size = legend_text_size),
              legend.spacing.x = unit(1.5, "cm"))
    )
    this_full_plot <- plot_grid(this_four_panel_plot, this_legend, ncol = 1, nrow = 2,
                                rel_heights = c(1, .1))
    plot_name <- paste0(
      "supplement_highdim_", vim, "_", this_corr
    )
    ggsave(filename = here("plots", paste0(plot_name, ".png")),
           plot = this_full_plot, device = "png", width = fig_width, height = fig_height, dpi = 300)
    ggsave(filename = here("plots", paste0(plot_name, ".eps")),
           plot = this_full_plot, device = "eps", width = fig_width, height = fig_height, dpi = 300)
  }
}
