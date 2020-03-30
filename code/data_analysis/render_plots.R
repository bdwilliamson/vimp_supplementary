#!/usr/local/bin/Rscript

## plot variable importance

## --------------------------------------------------------------------
## set up
## --------------------------------------------------------------------
library("dplyr")
library("tidyr")
library("vimp")
library("cowplot")
library("magick")
library("argparse")

source("../code/utils.R")
source("../code/get_variable_groups.R")
plots_dir <- "../code/plots/"

parser <- ArgumentParser()
parser$add_argument("--outcome", default = "cens", help = "outcome string")
parser$add_argument("--reduce-groups", action = "store_true", help = "should we run only for CD4 binding sites, or for all groups?")
parser$add_argument("--reduce-covs", action = "store_true", help = "should we run only for only a subset of the covariates?")
parser$add_argument("--reduce-library", action = "store_true", help = "should we run with a small learner library, or with all candidate learners?")
parser$add_argument("--threshold", type = "double", default = 0.05, help = "p-value threshold for 'significance'")
parser$add_argument("--point-size", type = "double", default = 3, help = "point size")
parser$add_argument("--axis-text-size", type = "double", default = 18, help = "axis text size")
parser$add_argument("--plot-text-size", type = "double", default = 6, help = "plot text size")
parser$add_argument("--main-text-size", type = "double", default = 20, help = "main text size")
args <- parser$parse_args()
print(args)

## read in the data
dat <- readRDS("../code/data/analysis_data.rds")

## if reduce_library, run with small lib
source("../code/super_learner_libraries.R")
if (args$reduce_library) {
    learner_lib <- library_reduced
} else {
    learner_lib <- library_full
}

## get names of predictors
all_var_groups <- get_variable_groups(dat)
## get variable groups
all_geog_vars <- all_var_groups$geog
## get all predictor names
pred_names <- unique(do.call(c, all_var_groups))
# if reduce_groups, only run on the cd4 binding site
if (args$reduce_groups) {
    var_groups <- all_var_groups[1]
} else {
    var_groups <- all_var_groups
}
# get names of outcomes
all_outcome_names <- c("ic50.censored", "binding.dichotomous.sens.resis")
outcome_names <- all_outcome_names[grepl(args$outcome, all_outcome_names)]
dat$ic50.censored <- as.numeric(dat$ic50.censored)
## specify the individual variables
if (args$reduce_covs) {
    num_covs <- 10
    var_inds <- pred_names[!grepl("geog", pred_names)][1:(num_covs - length(all_geog_vars))]
} else {
    num_covs <- length(pred_names) - length(all_geog_vars)
    var_inds <- pred_names[!grepl("geog", pred_names)][1:num_covs]
}

## --------------------------------------------------------------------
## load in the variable importance objects
## --------------------------------------------------------------------
num_pop_import <- 20
# x_lab_rsq <- expression(paste("Difference in ", R^2, sep = ""))
# x_lab_auc <- expression(paste("Difference in ", AUC, sep = ""))
# x_lab_acc <- expression(paste("Difference in ", "Classification Accuracy", sep = ""))
x_lab_rsq <- x_lab_auc <- x_lab_acc <- "Estimated Variable Importance"
vimp_measures <- c("r_squared", "accuracy", "auc")
x_lim_cond <- c(0, 0.4)
x_lim_indi <- c(0, 0.6)
x_lim_marg <- c(0, 0.75)
x_lim_lst <- list(x_lim_cond, x_lim_marg, x_lim_indi)

## plot options
fig_width <- 14
fig_height <- 14

imp_nms <- list(var_groups, var_groups, var_inds)
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    for (j in 1:length(vimp_measures)) {
        this_vimp_measure <- vimp_measures[j]
        if (this_vimp_measure == "r_squared") {
            this_x_lab <- x_lab_rsq
            this_plot_title <- expression(bold(R^2))
            this_x_lim_marg <- c(0, 1.4)
            this_x_lim_cond <- c(0, 1)
            this_x_lim_indi <- c(0, 1)
            this_x_lim_lst <- list(this_x_lim_cond, this_x_lim_marg, this_x_lim_indi)
            this_axis_font_size <- 14
        } else if (this_vimp_measure == "accuracy") {
            this_x_lab <- x_lab_acc
            this_plot_title <- expression(bold(Accuracy))
            this_x_lim_lst <- x_lim_lst
            this_axis_font_size <- args$axis_text_size
        } else {
            this_x_lab <- x_lab_auc
            this_x_lim_cond <- c(0, 0.6)
            this_x_lim_marg <- c(0, 0.6)
            this_x_lim_indi <- c(0, 0.6)
            this_plot_title <- expression(bold(AUC))
            this_x_lim_lst <- list(this_x_lim_cond, this_x_lim_marg, this_x_lim_indi)
            this_axis_font_size <- args$axis_text_size
        }
        ## importance results
        eval(parse(text = paste0(this_outcome_name, "_", this_vimp_measure, "_cv_vimp_lst <- readRDS(file = '../code/slfits/", this_outcome_name, "_", this_vimp_measure, "_cv_vimp.rds')")))
        ## make plots
        eval(parse(text = paste0("current_cv_vimp_lst <- ", this_outcome_name, "_", this_vimp_measure, "_cv_vimp_lst")))
        # vimp_plot_titles <- paste0(vimp_plot_name(this_outcome_name), ": ", names(current_cv_vimp_lst))
        plot_title_expr <- vimp_plot_name_expr(this_outcome_name)
        vimp_plot_titles <- rep(list(eval(bquote(expression(.(plot_title_expr[[1]]) ~ .(this_plot_title[[1]]))))), 3)
        current_vimp_plots <- mapply(function(x, y, z) plot_one_vimp(x, title = y, x_lab = this_x_lab, x_lim = z, cv = TRUE, num_plot = num_pop_import, threshold = args$threshold, point_size = args$point_size, main_font_size = args$main_text_size, axis_font_size = this_axis_font_size, plot_font_size = args$plot_text_size), current_cv_vimp_lst, vimp_plot_titles, this_x_lim_lst, SIMPLIFY = FALSE)
        eval(parse(text = paste0(this_outcome_name, "_", this_vimp_measure, "_vimp_plots <- current_vimp_plots")))
        ## conditional, marginal, individual plots
        ggsave(filename = paste0(plots_dir, this_outcome_name, '_', this_vimp_measure, '_conditional.png'), plot = current_vimp_plots$conditional, width = fig_width, height = fig_height, units = 'cm')
        ggsave(filename = paste0(plots_dir, this_outcome_name, '_', this_vimp_measure, '_marginal.png'), plot = current_vimp_plots$marginal, width = fig_width, height = fig_height, units = 'cm')
        ggsave(filename = paste0(plots_dir, this_outcome_name, '_', this_vimp_measure, '_individual.png'), plot = current_vimp_plots$individual, width = fig_width, height = fig_height, units = 'cm')
        
        ## a second round of plots with simpler plot title
        plot_title_simple <- vimp_plot_name_expr_simple(this_outcome_name)
        vimp_plot_titles_simple <- rep(list(eval(bquote(expression(.(plot_title_simple[[1]]) ~ .(this_plot_title[[1]]))))), 3)
        current_vimp_plots_simple <- mapply(function(x, y, z) plot_one_vimp(x, title = y, x_lab = this_x_lab, x_lim = z, cv = TRUE, num_plot = num_pop_import, threshold = args$threshold, point_size = args$point_size, main_font_size = args$main_text_size, axis_font_size = this_axis_font_size, plot_font_size = args$plot_text_size), current_cv_vimp_lst, vimp_plot_titles_simple, this_x_lim_lst, SIMPLIFY = FALSE)
        eval(parse(text = paste0(this_outcome_name, "_", this_vimp_measure, "_vimp_plots_simple <- current_vimp_plots_simple")))
        ## conditional, marginal, individual plots
        ggsave(filename = paste0(plots_dir, this_outcome_name, '_', this_vimp_measure, '_conditional_simple.png'), plot = current_vimp_plots_simple$conditional, width = fig_width, height = fig_height, units = 'cm')
        ggsave(filename = paste0(plots_dir, this_outcome_name, '_', this_vimp_measure, '_marginal_simple.png'), plot = current_vimp_plots_simple$marginal, width = fig_width, height = fig_height, units = 'cm')
        ggsave(filename = paste0(plots_dir, this_outcome_name, '_', this_vimp_measure, '_individual_simple.png'), plot = current_vimp_plots_simple$individual, width = fig_width, height = fig_height, units = 'cm')
        
    }
}

## --------------------------------------------------------------------
## make the plots
## --------------------------------------------------------------------
width_mult <- 1
outer_fig_width <- 50

## create feature group definitions
top_row <- 0.875
spacing <- 0.05
num_word_space <- 0.1
title_font_size <- 18
main_font_size <- 16
feature_group_defs <- ggdraw() +
  draw_label("Feature group definitions", size = title_font_size, x = 0.5, y = 0.95) +
  geom_hline(yintercept = 0.9) +
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
  draw_label("14", size = main_font_size, y = top_row - 15*spacing, x = 0, hjust = 0) +
  draw_label("VRC01 binding footprint", size = main_font_size, y = top_row, x = num_word_space, hjust = 0.0) +
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
  draw_label("Steric bulk at critical locations", size = main_font_size, y = top_row - 14*spacing, x = num_word_space, hjust = 0) +
  draw_label("Geographic confounders", size = main_font_size, y = top_row - 15*spacing, x = num_word_space, hjust = 0)

ggsave(filename = paste0(plots_dir, "vim_grp_definitions.png"),
     plot = feature_group_defs,
     width = width_mult*outer_fig_width/3, height = outer_fig_width/3, units = "cm")


## conditional, marginal, individual plots of accuracy, auc, and feature groups definitions
for (i in 1:length(outcome_names)) {
    this_outcome_name <- outcome_names[i]
    ## first, the one for the paper
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_acc_auc_conditional.png'), plot = plot_grid(", this_outcome_name, "_accuracy_vimp_plots$conditional, ", this_outcome_name, "_auc_vimp_plots$conditional, feature_group_defs, labels = 'AUTO', nrow = 1), width = width_mult*outer_fig_width/1.25, height = outer_fig_width/3, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_acc_auc_marginal.png'), plot = plot_grid(", this_outcome_name, "_accuracy_vimp_plots$marginal, ", this_outcome_name, "_auc_vimp_plots$marginal, feature_group_defs, labels = 'AUTO', nrow = 1), width = width_mult*outer_fig_width/1.25, height = outer_fig_width/3, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_acc_auc_individual.png'), plot = plot_grid(", this_outcome_name, "_accuracy_vimp_plots$individual, ", this_outcome_name, "_auc_vimp_plots$individual, labels = 'AUTO', nrow = 1), width = width_mult*outer_fig_width/1.25, height = outer_fig_width/3, units = 'cm')")))

    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_r2_conditional.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$conditional, feature_group_defs, labels = 'AUTO', nrow = 1), width = width_mult*outer_fig_width/2, height = outer_fig_width/3, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_r2_marginal.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$marginal,  feature_group_defs, labels = 'AUTO', nrow = 1), width = width_mult*outer_fig_width/2, height = outer_fig_width/3, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_r2_individual.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$individual, labels = 'AUTO', nrow = 1), width = width_mult*outer_fig_width/2, height = outer_fig_width/3, units = 'cm')")))
}

## make plots for presentations
for (i in 1:length(outcome_names)) {
    ## second, the one for presentations
    this_outcome_name <- outcome_names[i]
    ## r-squared only
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_conditional.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$conditional, plot.new(), plot.new(), feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_marginal.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$marginal, plot.new(), plot.new(), feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_individual.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$individual, plot.new(), plot.new(), labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
    ## with accuracy and auc as well
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_acc_auc_conditional.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$conditional, ", this_outcome_name, "_accuracy_vimp_plots$conditional, ", this_outcome_name, "_auc_vimp_plots$conditional, feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_acc_auc_marginal.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots$marginal, ", this_outcome_name, "_accuracy_vimp_plots$marginal, ", this_outcome_name, "_auc_vimp_plots$marginal, feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
    eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_acc_auc_individual.png'), plot = plot_grid(",this_outcome_name, "_r_squared_vimp_plots$individual, ", this_outcome_name, "_accuracy_vimp_plots$individual, ", this_outcome_name, "_auc_vimp_plots$individual, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
}

## make plots for presentations with simplified outcome "resistant" means ic50 was censored (i.e., concentration too high to be measured)
for (i in 1:length(outcome_names)) {
  ## second, the one for presentations
  this_outcome_name <- outcome_names[i]
  ## r-squared only
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_conditional_simple.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots_simple$conditional, plot.new(), plot.new(), feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_marginal_simple.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots_simple$marginal, plot.new(), plot.new(), feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_individual_simple.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots_simple$individual, plot.new(), plot.new(), labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
  ## with accuracy and auc as well
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_acc_auc_conditional_simple.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots_simple$conditional, ", this_outcome_name, "_accuracy_vimp_plots_simple$conditional, ", this_outcome_name, "_auc_vimp_plots_simple$conditional, feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_acc_auc_marginal_simple.png'), plot = plot_grid(", this_outcome_name, "_r_squared_vimp_plots_simple$marginal, ", this_outcome_name, "_accuracy_vimp_plots_simple$marginal, ", this_outcome_name, "_auc_vimp_plots_simple$marginal, feature_group_defs, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'vim_', this_outcome_name, '_pres_r2_acc_auc_individual_simple.png'), plot = plot_grid(",this_outcome_name, "_r_squared_vimp_plots_simple$individual, ", this_outcome_name, "_accuracy_vimp_plots_simple$individual, ", this_outcome_name, "_auc_vimp_plots_simple$individual, labels = 'AUTO', label_size = 25), width = width_mult*outer_fig_width/1.5, height = outer_fig_width/2, units = 'cm')")))
}

## --------------------------------------------------------------------------
## Super Learner predictor descriptions, for manuscript
## --------------------------------------------------------------------------
## table of the learner library and weights
sl_descrip_table <- sapply(learner_lib, function(x){
	eval(parse(text = paste0("descr_", x)))
})
sl_descrip_tbl <- tibble::tibble(function_name = names(sl_descrip_table), description = sl_descrip_table)
knitr::kable(sl_descrip_tbl, format = "latex", col.names = c("Function name", "Description"), caption = "Library of canididate learners for the Super Learner with descriptions.", label = "sl_lib") %>%
    cat(., file = "../code/plots/sl_lib_table.tex")
# load super learner fits
weight_list <- vector(mode = "list", length = length(outcome_names))
ct <- 0
for(o in outcome_names){
	ct <- ct + 1
	fit_name <- paste0("../code/slfits/cvfit_", o, ".rds")
	fit <- readRDS(fit_name)
	class(fit) <- "myCV.SuperLearner"
	weight_list[[ct]] <- fit$coef
}
weight_table <- data.frame(Reduce(cbind, weight_list))
knitr::kable(as.data.frame(weight_table), format = "latex", digits = 2, caption = "Table of Super Learner weights for each candidate learner and cross-validation fold.", label = "sl_weights") %>%
  cat(., file = "../code/plots/sl_weight_table.tex")

## plot of cv performance of the super learner in terms of AUC
for (o in outcome_names) {
  fit_name <- paste0("../code/slfits/cvfit_", o, ".rds")
  fit <- readRDS(fit_name)
  class(fit) <- "myCV.SuperLearner"
  main_title <- ifelse(o == "ic50.censored", expression(bold(paste(IC[50], " Censored"))), expression(bold(paste("Sensitive/Resistant Only"))))
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'sl_perf_', o, '.png'), plot = plot(fit, main_title = main_title, xlim1 = 0.4, xlim2 = 1, text_size = 9, Rsquared = FALSE, method = 'method.AUC'), width = fig_width, height = fig_height, units = 'cm')")))
  eval(parse(text = paste0("ggsave(filename = paste0(plots_dir, 'sl_roc_', o, '.png'), plot = plot_roc_curves(fit), width = fig_width * 1.5, height = fig_height, units = 'cm')")))
}
