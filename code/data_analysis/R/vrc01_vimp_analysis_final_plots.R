## produce paper-ready plots based on the individual plots from separate vimp measures

## ---------------------------------
## setup
## ---------------------------------
## load required functions and packages
library("magick")
library("cowplot")

## ---------------------------------
## read in the plots
## ---------------------------------
## setup
plots_dir <- "../code/plots/"
fig_width <- 50
outcomes <- c("cens", "sens.resis")
measures <- c("accuracy", "auc", "deviance")
outcome_measure_df <- expand.grid(measures, outcomes)
outcome_lst <- as.list(outcome_measure_df[, 2])
measure_lst <- as.list(outcome_measure_df[, 1])

## read in all of the plots
all_plots <- mapply(function(x, y) image_read(paste0(plots_dir, "vim_grp_combined_outcome_", x, "_standard_level_0.95_", y, ".png")), outcome_lst, measure_lst)

## split off into separate lists for each measure
cens_plots <- all_plots[1:3]
sens_resis_plots <- all_plots[4:6]

## ------------------------------------------------------
## create new plots: one for censored, one for sens/resis
## ------------------------------------------------------
## create feature group definitions
top_row <- 0.875
spacing <- 0.06
num_word_space <- 0.1
title_font_size <- 13
main_font_size <- 11
feature_group_defs <- ggdraw() +
  draw_label("Feature group definitions", size = title_font_size, x = 0.5, y = 0.95) +
  geom_hline(yintercept = 0.92) +
  draw_label("1", size = main_font_size, y = top_row, x = 0, hjust = 0) +
  draw_label("2", size = main_font_size, y = top_row - spacing, x = 0, hjust = 0) +
  draw_label("3", size = main_font_size, y = top_row - 2*spacing, x = 0, hjust = 0) +
  draw_label("4", size = main_font_size, y = top_row - 3*spacing, x = 0, hjust = 0) +
  draw_label("5", size = main_font_size, y = top_row - 4*spacing, x = 0, hjust = 0) +
  draw_label("", size = main_font_size, y = top_row - 5*spacing, x = 0, hjust = 0) +
  draw_label("6", size = main_font_size, y = top_row - 6*spacing, x = 0, hjust = 0) +
  draw_label("", size = main_font_size, y = top_row - 7*spacing, x = 0, hjust = 0) +
  draw_label("7", size = main_font_size, y = top_row - 8*spacing, x = 0, hjust = 0) +
  draw_label("8", size = main_font_size, y = top_row - 9*spacing, x = 0, hjust = 0) +
  draw_label("9", size = main_font_size, y = top_row - 10*spacing, x = 0, hjust = 0) +
  draw_label("10", size = main_font_size, y = top_row - 11*spacing, x = 0, hjust = 0) +
  draw_label("11", size = main_font_size, y = top_row - 12*spacing, x = 0, hjust = 0) +
  draw_label("12", size = main_font_size, y = top_row - 13*spacing, x = 0, hjust = 0) +
  draw_label("13", size = main_font_size, y = top_row - 14*spacing, x = 0, hjust = 0) +
  draw_label("VRC01 binding footprint", size = main_font_size, y = 0.875, x = num_word_space, hjust = 0.0) +
  draw_label("CD4 binding sites", size = main_font_size, y = top_row - spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites with sufficient exposed surface area", size = main_font_size, y = top_row - 2*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites identified as important for glycosylation", size = main_font_size, y = top_row - 3*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites with residues that covary with the ", size = main_font_size, y = top_row - 4*spacing, x = num_word_space, hjust = 0) +
  draw_label("VRC01 binding footprint", size = main_font_size, y = top_row - 5*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites associated with VRC01-specific ", size = main_font_size, y = top_row - 6*spacing, x = num_word_space, hjust = 0) +
  draw_label("potential N-linked glycosylation (PNGS) effects", size = main_font_size, y = top_row - 7*spacing, x = num_word_space, hjust = 0) +
  draw_label("gp41 sites important for VRC01 binding", size = main_font_size, y = top_row - 8*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites for indicating N-linked glycosylation", size = main_font_size, y = top_row - 9*spacing, x = num_word_space, hjust = 0) +
  draw_label("Majority virus subtypes", size = main_font_size, y = top_row - 10*spacing, x = num_word_space, hjust = 0) +
  draw_label("Region-specific counts of PNGS", size = main_font_size, y = top_row - 11*spacing, x = num_word_space, hjust = 0) +
  draw_label("Viral geometry", size = main_font_size, y = top_row - 12*spacing, x = num_word_space, hjust = 0) +
  draw_label("Cysteine counts", size = main_font_size, y = top_row - 13*spacing, x = num_word_space, hjust = 0) +
  draw_label("Steric bulk at critical locations", size = main_font_size, y = top_row - 14*spacing, x = num_word_space, hjust = 0)
## ----------------
## create cens plot
## ----------------
cens_1 <- ggdraw() + draw_image(cens_plots[[1]])
cens_2 <- ggdraw() + draw_image(cens_plots[[2]])
cens_3 <- ggdraw() + draw_image(cens_plots[[3]])
width_mult <- 1

ggsave(filename = paste0(plots_dir,         "vim_grp_outcome_cens_standard_level_0.95.png"),
plot = plot_grid(cens_1, cens_2, cens_3, feature_group_defs, labels = "AUTO"),
width = width_mult*fig_width, height = fig_width, units = "cm")

ggsave(filename = paste0(plots_dir,         "vim_grp_outcome_cens_standard_level_0.95.tiff"),
plot = plot_grid(cens_1, cens_2, cens_3, feature_group_defs, labels = "AUTO"),
width = width_mult*fig_width, height = fig_width, units = "cm")

## ----------------
## create sens/resis plot
## ----------------
sens_resis_1 <- ggdraw() + draw_image(sens_resis_plots[[1]])
sens_resis_2 <- ggdraw() + draw_image(sens_resis_plots[[2]])
sens_resis_3 <- ggdraw() + draw_image(sens_resis_plots[[3]])

ggsave(filename = paste0(plots_dir, "vim_grp_outcome_sens.resis_standard_level_0.95.png"), plot = plot_grid(sens_resis_1, sens_resis_2, sens_resis_3, feature_group_defs, labels = "AUTO"),
width = width_mult*fig_width, height = fig_width, units = "cm")

ggsave(filename = paste0(plots_dir, "vim_grp_outcome_sens.resis_standard_level_0.95.tiff"), plot = plot_grid(sens_resis_1, sens_resis_2, sens_resis_3, feature_group_defs, labels = "AUTO"),
width = width_mult*fig_width, height = fig_width, units = "cm")
