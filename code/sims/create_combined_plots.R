#!/usr/local/bin/Rscript

# create second-level plots for supplement

library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("here")

#### set up objects ####
all_vims <- c("accuracy", "auc")
alt_js <- c(1, 2)
null_js <- c(2, 3)
all_alt <- apply(expand.grid(vim = all_vims, j = alt_js), 1, function(x) paste0(x[1], "_", x[2]))
all_null <- apply(expand.grid(vim = all_vims, j = null_js), 1, function(x) paste0(x[1], "_", x[2]))
alt_names <- paste0("supplement_alternative_", all_alt, "_p_2_boot_0_nocf.rds")
null_names <- paste0("supplement_null_", all_null, "_p_4_boot_0_nocf.rds")

# read in the plots
alt_plots_lst <- lapply(
    as.list(here("plots", alt_names)), readRDS
)
null_plots_lst <- lapply(
    as.list(here("plots", null_names)), readRDS
)
# note that the legends are all the same
legend <- readRDS(here("plots", "supplement_alternative_accuracy_2_p_2_boot_0_nocf_legend.rds"))

# useful objects
point_size <- 1.5
title_text_size <- 12
legend_text_size <- 10
axis_text_size <- 10
fig_width <- 9
fig_height <- 4.25 * 1.5
dodge_width <- 0.875
color_begin <- 0
color_end <- 0.5

#### create a 4x4 plot for the alternative: coverage ####
# first, remove the legend from all but the first one
alt_combined_plot <- plot_grid(
    plotlist = alt_plots_lst, labels = "auto"
)
alt_final_plot <- plot_grid(alt_combined_plot, legend, ncol = 1, nrow = 2,
                            rel_heights = c(1, .1))
ggsave(filename = here("plots", "supplement_alternative_combined_cover.png"),
       alt_final_plot, width = fig_width, height = fig_height)
ggsave(filename = here("plots", "supplement_alternative_combined_cover.eps"),
       alt_final_plot, width = fig_width, height = fig_height)

#### create a 4x4 plot for the null: coverage and type 1 error ####
# first, remove the legend from all but the first one
null_combined_plot <- plot_grid(
    plotlist = null_plots_lst, labels = "auto"
)
null_final_plot <- plot_grid(null_combined_plot, legend, ncol = 1, nrow = 2,
                             rel_heights = c(1, .1))
ggsave(filename = here("plots", "supplement_null_combined_cover.png"),
       null_combined_plot, width = fig_width * 1.75, height = fig_height)
ggsave(filename = here("plots", "supplement_null_combined_cover.eps"),
       null_combined_plot, width = fig_width * 1.75, height = fig_height)
