#!/usr/local/bin/Rscript

# plot variable importance

#### set up ####
library("dplyr")
library("tidyr")
library("vimp")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("argparse")
library("here")
library("kableExtra")

source(here("code", "utils.R"))
source(here("code", "get_variable_groups.R"))
plots_dir <- "plots/all_vrc01_plots"
if (!dir.exists(here(plots_dir))) {
    dir.create(here(plots_dir), recursive = TRUE)
}

parser <- ArgumentParser()
parser$add_argument("--outcome", default = c("cens", "sens50"),
                    help = paste0("outcome string: ",
                                  "cens for right-censored IC-50,",
                                  "sens50 for IC-50 < 1,",
                                  "sens80 for IC-80 < 1."),
                    nargs = "+"
)
parser$add_argument("--ind-type", default = "sitewise",
                    help = "should we do site-wise or residue-wise individual importance?")
parser$add_argument("--reduce-groups", action = "store_true",
                    help = "should we run only for CD4 binding sites, or for all groups?")
parser$add_argument("--reduce-covs", action = "store_true",
                    help = "should we run only for only a subset of the covariates?")
parser$add_argument("--reduce-library", action = "store_true",
                    help = "should we run with a small learner library, or with all candidate learners?")
parser$add_argument("--threshold", type = "double", default = 0.05 / 13,
                    help = "p-value threshold for 'significance'")
parser$add_argument("--point-size", type = "double", default = 2, help = "point size")
parser$add_argument("--axis-text-size", type = "double", default = 10,
                    help = "axis text size")
parser$add_argument("--plot-text-size", type = "double", default = 4,
                    help = "plot text size")
parser$add_argument("--main-text-size", type = "double", default = 12,
                    help = "main text size")
args <- parser$parse_args()
print(args)

# read in the data
dat <- readRDS(here("code", "data", "analysis_data.rds"))

# if reduce_library, run with small lib
source(here("code", "super_learner_libraries.R"))
if (args$reduce_library) {
    learner_lib <- library_reduced
} else {
    learner_lib <- library_full
}

# specify the individual variables and groups
# get names of predictors
all_var_groups <- get_variable_groups(dat)
# get variable groups
all_geog_vars <- all_var_groups$geog
# get all predictor names
pred_names <- unique(do.call(c, all_var_groups))
if (args$reduce_covs) {
    num_covs <- 10
} else {
    num_covs <- length(pred_names)
}
# get all individual predictors
all_var_inds <- get_individual_features(pred_names[!grepl("geog", pred_names)][1:num_covs],
                                        args$ind_type)

# if reduce_groups, only run on the cd4 binding site
if (args$reduce_groups) {
    var_groups <- all_var_groups[1]
} else {
    var_groups <- all_var_groups
}
# the only possible outcomes for these analyses are: IC50 censored ("ic50.censored"),
# sensitive/resistant ("binding.dichotomous.sens.resis"), IC50 < 1, and IC80 < 1
initial_analysis_dataset <- dat[, c("ic50.censored", "binding.dichotomous.sens.resis",
                                    "sens50", "sens80", pred_names)]

#### read in variable importance objects and make individual plots ####
num_pop_import <- 20
x_lab_rsq <- x_lab_auc <- x_lab_acc <- "Estimated Variable Importance"
vimp_measures <- c("r_squared", "accuracy", "auc")
x_pad <- 0.35

# plot options
fig_width <- 9
fig_height <- 4.25
filetypes <- c("png", "tiff", "eps")

imp_nms <- list(var_groups, var_groups, all_var_inds)
for (i in seq_len(length(args$outcome))) {
    this_outcome_name <- get_outcome(args$outcome[i])
    for (j in 1:length(vimp_measures)) {
        this_vimp_measure <- vimp_measures[j]
        if (this_vimp_measure == "r_squared") {
            this_x_lab <- x_lab_rsq
            this_plot_title <- expression(bold(R^2))
        } else if (this_vimp_measure == "accuracy") {
            this_x_lab <- x_lab_acc
            this_plot_title <- expression(bold(Accuracy))
        } else {
            this_x_lab <- x_lab_auc
            this_plot_title <- expression(bold(AUC))
        }
        this_axis_font_size <- args$axis_text_size
        # importance results
        eval(parse(
            text = paste0(
                this_outcome_name, "_", this_vimp_measure,
                "_cv_vimp_lst <- readRDS(file = here('results', '", this_outcome_name,
                "_", this_vimp_measure, "_cv_vimp.rds'))"
            )
        ))
        print(paste0("Outcome: ", this_outcome_name, "; ",
                     "VIM: ", this_vimp_measure))
        print("Conditional p-values:")
        eval(parse(text = paste0("print(", this_outcome_name, "_", this_vimp_measure,
                                 "_cv_vimp_lst$conditional$mat %>%
                                 dplyr::filter(p_value < 0.05) %>%
                                 dplyr::select(s, p_value))")))
        print("Marginal p-values:")
        eval(parse(text = paste0("print(", this_outcome_name, "_", this_vimp_measure,
                                 "_cv_vimp_lst$marginal$mat %>%
                                 dplyr::filter(p_value < 0.05) %>%
                                 dplyr::select(s, p_value))")))
        # make plots
        eval(parse(
            text = paste0("current_cv_vimp_lst <- ", this_outcome_name, "_",
                          this_vimp_measure, "_cv_vimp_lst")
        ))
        this_x_lim_lst <- list(
            c(0, max(current_cv_vimp_lst$conditional$ci) + x_pad),
            c(0, max(current_cv_vimp_lst$marginal$ci) + x_pad),
            c(0, max(current_cv_vimp_lst$individual$ci) + x_pad)
        )
        plot_title_expr <- vimp_plot_name_expr(this_outcome_name)
        vimp_plot_titles <- rep(list(
            eval(bquote(expression(.(plot_title_expr[[1]]) ~ .(this_plot_title[[1]]))))
        ), 3)
        current_vimp_plots <- mapply(
            function(x, y, z, w) {
                plot_one_vimp(x, title = y, x_lab = this_x_lab, x_lim = z, cv = TRUE,
                              num_only = w,
                              num_plot = num_pop_import, threshold = args$threshold,
                              point_size = args$point_size,
                              main_font_size = args$main_text_size,
                              axis_font_size = this_axis_font_size,
                              plot_font_size = args$plot_text_size)
                },
            current_cv_vimp_lst, vimp_plot_titles, this_x_lim_lst,
            !grepl("individual", names(current_cv_vimp_lst)), SIMPLIFY = FALSE
        )
        eval(parse(
            text = paste0(this_outcome_name, "_", this_vimp_measure,
                          "_vimp_plots <- current_vimp_plots")
        ))
        # a second round of plots with simpler plot title
        plot_title_simple <- vimp_plot_name_expr_simple(this_outcome_name)
        vimp_plot_titles_simple <- rep(list(
            eval(bquote(expression(.(plot_title_simple[[1]]) ~ .(this_plot_title[[1]]))))
        ), 3)
        current_vimp_plots_simple <- mapply(
            function(x, y, z, w) {
                plot_one_vimp(x, title = y, x_lab = this_x_lab, x_lim = z, cv = TRUE,
                              num_only = w,
                              num_plot = num_pop_import, threshold = args$threshold,
                              point_size = args$point_size,
                              main_font_size = args$main_text_size,
                              axis_font_size = this_axis_font_size,
                              plot_font_size = args$plot_text_size)
            }, current_cv_vimp_lst, vimp_plot_titles_simple, this_x_lim_lst,
            !grepl("individual", names(current_cv_vimp_lst)), SIMPLIFY = FALSE)
        eval(parse(
            text = paste0(this_outcome_name, "_", this_vimp_measure,
                          "_vimp_plots_simple <- current_vimp_plots_simple")
        ))
        # conditional, marginal, individual plots
        filename_common <- paste0(this_outcome_name, "_", this_vimp_measure)
        for (filetype in filetypes) {
            ggsave(filename = here(plots_dir, paste0(filename_common, "_conditional.", filetype)),
                   plot = current_vimp_plots$conditional,
                   width = fig_width, height = fig_height, units = 'in', dpi = 300)
            ggsave(filename = here(plots_dir, paste0(filename_common, "_marginal.", filetype)),
                   plot = current_vimp_plots$marginal,
                   width = fig_width, height = fig_height, units = 'in', dpi = 300)
            ggsave(filename = here(plots_dir, paste0(filename_common, "_individual.", filetype)),
                   plot = current_vimp_plots$individual,
                   width = fig_width, height = fig_height, units = 'in', dpi = 300)

            ggsave(filename = here(plots_dir, paste0(filename_common, "_conditional_simple.", filetype)),
                   plot = current_vimp_plots_simple$conditional,
                   width = fig_width, height = fig_height, units = 'in', dpi = 300)
            ggsave(filename = here(plots_dir, paste0(filename_common, "_marginal_simple.", filetype)),
                   plot = current_vimp_plots_simple$marginal,
                   width = fig_width, height = fig_height, units = 'in', dpi = 300)
            ggsave(filename = here(plots_dir, paste0(filename_common, "_individual_simple.", filetype)),
                   plot = current_vimp_plots_simple$individual,
                   width = fig_width, height = fig_height, units = 'in', dpi = 300)
        }
    }
}

#### make a plot with the feature group definitions ####
width_mult <- 1
outer_fig_width <- 9
outer_fig_height <- fig_height

# create feature group definitions
top_row <- 0.875
spacing <- 0.05
num_word_space <- 0.1
feature_grp_text_size <- args$main_text_size * 0.75
feature_group_defs <- ggdraw() +
  draw_label("Feature group definitions", size = args$main_text_size, x = 0.5, y = 0.95) +
  geom_hline(yintercept = 0.9) +
  draw_label("1", size = feature_grp_text_size, y = top_row, x = 0, hjust = 0) +
  draw_label("2", size = feature_grp_text_size, y = top_row - spacing, x = 0, hjust = 0) +
  draw_label("3", size = feature_grp_text_size, y = top_row - 2*spacing, x = 0, hjust = 0) +
  draw_label("4", size = feature_grp_text_size, y = top_row - 3*spacing, x = 0, hjust = 0) +
  draw_label("5", size = feature_grp_text_size, y = top_row - 4*spacing, x = 0, hjust = 0) +
  draw_label("", size = feature_grp_text_size, y = top_row - 5*spacing, x = 0, hjust = 0) +
  draw_label("6", size = feature_grp_text_size, y = top_row - 6*spacing, x = 0, hjust = 0) +
  draw_label("", size = feature_grp_text_size, y = top_row - 7*spacing, x = 0, hjust = 0) +
  draw_label("7", size = feature_grp_text_size, y = top_row - 8*spacing, x = 0, hjust = 0) +
  draw_label("8", size = feature_grp_text_size, y = top_row - 9*spacing, x = 0, hjust = 0) +
  draw_label("9", size = feature_grp_text_size, y = top_row - 10*spacing, x = 0, hjust = 0) +
  draw_label("10", size = feature_grp_text_size, y = top_row - 11*spacing, x = 0, hjust = 0) +
  draw_label("11", size = feature_grp_text_size, y = top_row - 12*spacing, x = 0, hjust = 0) +
  draw_label("12", size = feature_grp_text_size, y = top_row - 13*spacing, x = 0, hjust = 0) +
  draw_label("13", size = feature_grp_text_size, y = top_row - 14*spacing, x = 0, hjust = 0) +
  draw_label("14", size = feature_grp_text_size, y = top_row - 15*spacing, x = 0, hjust = 0) +
  draw_label("VRC01 binding footprint", size = feature_grp_text_size, y = top_row, x = num_word_space, hjust = 0.0) +
  draw_label("CD4 binding sites", size = feature_grp_text_size, y = top_row - spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites with sufficient exposed surface area", size = feature_grp_text_size, y = top_row - 2*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites identified as important for glycosylation", size = feature_grp_text_size, y = top_row - 3*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites with residues that covary with the ", size = feature_grp_text_size, y = top_row - 4*spacing, x = num_word_space, hjust = 0) +
  draw_label("VRC01 binding footprint", size = feature_grp_text_size, y = top_row - 5*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites associated with VRC01-specific ", size = feature_grp_text_size, y = top_row - 6*spacing, x = num_word_space, hjust = 0) +
  draw_label("potential N-linked glycosylation (PNGS) effects", size = feature_grp_text_size, y = top_row - 7*spacing, x = num_word_space, hjust = 0) +
  draw_label("gp41 sites important for VRC01 binding", size = feature_grp_text_size, y = top_row - 8*spacing, x = num_word_space, hjust = 0) +
  draw_label("Sites for indicating N-linked glycosylation", size = feature_grp_text_size, y = top_row - 9*spacing, x = num_word_space, hjust = 0) +
  draw_label("Majority virus subtypes", size = feature_grp_text_size, y = top_row - 10*spacing, x = num_word_space, hjust = 0) +
  draw_label("Region-specific counts of PNGS", size = feature_grp_text_size, y = top_row - 11*spacing, x = num_word_space, hjust = 0) +
  draw_label("Viral geometry", size = feature_grp_text_size, y = top_row - 12*spacing, x = num_word_space, hjust = 0) +
  draw_label("Cysteine counts", size = feature_grp_text_size, y = top_row - 13*spacing, x = num_word_space, hjust = 0) +
  draw_label("Steric bulk at critical locations", size = feature_grp_text_size, y = top_row - 14*spacing, x = num_word_space, hjust = 0) +
  draw_label("Geographic confounders", size = feature_grp_text_size, y = top_row - 15*spacing, x = num_word_space, hjust = 0)

for (filetype in filetypes) {
    if (length(args$outcome) == 1) {
        common_filename <- paste0("vimp_grp_definitions_", args$outcome, ".")
    } else {
        common_filename <- "vim_grp_definitions."
    }
    ggsave(filename = here(plots_dir, paste0(common_filename, filetype)),
           plot = feature_group_defs,
           width = fig_width, height = fig_height, units = "in", dpi = 300)
}

# conditional, marginal, individual plots of accuracy, auc, and feature groups definitions
ms_plots_dir <- "plots"
supp_r2_width_mult <- 3 / 4
for (i in seq_len(length(args$outcome))) {
    this_outcome_name <- get_outcome(args$outcome[i])
    main_filename <- paste0("main_vim_", this_outcome_name)
    supplemental_filename <- paste0("supplement_vim_", this_outcome_name)
    # first, the one for the paper
    for (filetype in filetypes) {
        # for the main manuscript: accuracy, auc, feature groups
        eval(parse(
            text = paste0("ggsave(filename = here(ms_plots_dir, paste0(main_filename,
                                                  '_conditional.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                                                   "_accuracy_vimp_plots$conditional, ",
                                                   this_outcome_name,
                                                   "_auc_vimp_plots$conditional,
                                                   feature_group_defs,
                                                   labels = 'AUTO', nrow = 1),
                                  width = outer_fig_width,
                                  height = outer_fig_height,
                                  units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(ms_plots_dir, paste0(main_filename,
                                                  '_marginal.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                                                   "_accuracy_vimp_plots$marginal, ",
                                                   this_outcome_name,
                                                   "_auc_vimp_plots$marginal,
                                                   feature_group_defs,
                                                   labels = 'AUTO', nrow = 1),
                                  width = outer_fig_width,
                                  height = outer_fig_height,
                                  units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(ms_plots_dir, paste0(main_filename,
                                                  '_individual.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                                                   "_accuracy_vimp_plots$individual, ",
                                                   this_outcome_name,
                                                   "_auc_vimp_plots$individual,
                                                   labels = 'AUTO', nrow = 1),
                                  width = outer_fig_width,
                                  height = outer_fig_height,
                                  units = 'in', dpi = 300)")
        ))
        # for the supplement: R-squared
        eval(parse(
            text = paste0("ggsave(filename = here(ms_plots_dir, paste0(supplemental_filename,
                                                  '_conditional.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                                                   "_r_squared_vimp_plots$conditional,
                                                   feature_group_defs,
                                                   labels = 'AUTO', nrow = 1),
                                  width = outer_fig_width * supp_r2_width_mult,
                                  height = outer_fig_height,
                                  units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(ms_plots_dir, paste0(supplemental_filename,
                                                  '_marginal.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                                                   "_r_squared_vimp_plots$marginal,
                                                   feature_group_defs,
                                                   labels = 'AUTO', nrow = 1),
                                  width = outer_fig_width * supp_r2_width_mult,
                                  height = outer_fig_height,
                                  units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(ms_plots_dir, paste0(supplemental_filename,
                                                  '_individual.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                                                   "_r_squared_vimp_plots$individual,
                                                   labels = 'AUTO', nrow = 1),
                                  width = outer_fig_width * supp_r2_width_mult,
                                  height = outer_fig_height,
                                  units = 'in', dpi = 300)")
        ))
    }
}

# make plots for presentations
pres_width_mult <- 1
pres_height_mult <- 2
for (i in seq_len(length(args$outcome))) {
    # second, the one for presentations
    this_outcome_name <- get_outcome(args$outcome[i])
    presentation_filename <- paste0("presentation_vim_", this_outcome_name)
    # r-squared only
    for (filetype in filetypes) {
        # r-squared only
        eval(parse(
            text = paste0("ggsave(filename = here(plots_dir,
                                paste0(presentation_filename, '_r2_conditional.', filetype)),
                              plot = plot_grid(", this_outcome_name,
                                 "_r_squared_vimp_plots_simple$conditional, plot.new(),
                                 plot.new(), feature_group_defs, labels = 'AUTO',
                                 label_size = 25),
                              width = outer_fig_width * pres_width_mult,
                              height = outer_fig_height * pres_height_mult,
                          units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(plots_dir,
                                paste0(presentation_filename, '_r2_marginal.', filetype)),
                              plot = plot_grid(", this_outcome_name,
                                 "_r_squared_vimp_plots_simple$marginal, plot.new(),
                                 plot.new(), feature_group_defs, labels = 'AUTO',
                                 label_size = 25),
                              width = outer_fig_width * pres_width_mult,
                              height = outer_fig_height * pres_height_mult,
                          units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(plots_dir,
                                paste0(presentation_filename, '_r2_individual.', filetype)),
                              plot = plot_grid(", this_outcome_name,
                                 "_r_squared_vimp_plots_simple$individual, plot.new(),
                                 plot.new(), labels = 'AUTO',
                                 label_size = 25),
                              width = outer_fig_width * pres_width_mult,
                              height = outer_fig_height * pres_height_mult,
                          units = 'in', dpi = 300)")
        ))
        # add in accuracy and AUC
        eval(parse(
            text = paste0("ggsave(filename = here(plots_dir,
                                    paste0(presentation_filename, '_all_conditional.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                                     "_r_squared_vimp_plots_simple$conditional, ",
                                     this_outcome_name, "_accuracy_vimp_plots_simple$conditional, ",
                                     this_outcome_name, "_auc_vimp_plots_simple$conditional,
                                     feature_group_defs, labels = 'AUTO', label_size = 25),
                                  width = outer_fig_width * pres_width_mult,
                              height = outer_fig_height * pres_height_mult,
                          units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(plots_dir,
                                    paste0(presentation_filename, '_all_marginal.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                          "_r_squared_vimp_plots_simple$marginal, ",
                          this_outcome_name, "_accuracy_vimp_plots_simple$marginal, ",
                          this_outcome_name, "_auc_vimp_plots_simple$marginal,
                                     feature_group_defs, labels = 'AUTO', label_size = 25),
                                  width = outer_fig_width * pres_width_mult,
                              height = outer_fig_height * pres_height_mult,
                          units = 'in', dpi = 300)")
        ))
        eval(parse(
            text = paste0("ggsave(filename = here(plots_dir,
                                    paste0(presentation_filename, '_all_individual.', filetype)),
                                  plot = plot_grid(", this_outcome_name,
                          "_r_squared_vimp_plots_simple$individual, ",
                          this_outcome_name, "_accuracy_vimp_plots_simple$individual, ",
                          this_outcome_name, "_auc_vimp_plots_simple$individual, labels = 'AUTO', label_size = 25),
                                  width = outer_fig_width * pres_width_mult,
                              height = outer_fig_height * pres_height_mult,
                          units = 'in', dpi = 300)")
        ))

    }
}

# --------------------------------------------------------------------------
# Super Learner predictor descriptions, for manuscript
# --------------------------------------------------------------------------
# table of the learner library and weights
sl_descrip_table <- sapply(learner_lib, function(x){
	eval(parse(text = paste0("descr_", x)))
})
sl_descrip_tbl <- tibble::tibble(function_name = names(sl_descrip_table),
                                 description = sl_descrip_table)
knitr::kable(sl_descrip_tbl, format = "latex", col.names = c("Function name", "Description"),
             caption = "Library of canididate learners for the Super Learner with descriptions. Unless otherwise specified, the default tuning parameters from \texttt{SuperLearner} are used. For boosted trees: a maximum of 1000 trees, a minimum of 10 observations per node, shrinkage parameter 0.1, logistic loss for binary outcomes and squared error loss for continuous outcomes. For random forests: 500 trees, minimum node size of 5 for continuous outcomes and 1 for binary outcomes, regression trees for continuous outcomes and probability trees for binary outcomes.", 
             label = "sl_lib") %>%
    kableExtra::kable_styling(latex_options = c("scale_down")) %>%
    cat(., file = here(ms_plots_dir, "supplement_sl_lib_table.tex"))
# load super learner fits
weight_list <- vector(mode = "list", length = length(args$outcome))
ct <- 0
for(o in args$outcome) {
	ct <- ct + 1
	if (o == "cens") {
	    o <- "ic50.censored"
	}
	fit_name <- paste0("cvfit_", o, ".rds")
	fit <- readRDS(here("code", "slfits", fit_name))
	class(fit) <- "myCV.SuperLearner"
	weight_list[[ct]] <- fit$coef
}
weight_table <- data.frame(Reduce(rbind, weight_list))
weight_table$fold <- rep(1:10, length(args$outcome))
colnames(weight_table) <- gsub("_All", "", gsub("SL.", "", colnames(weight_table)))
if (length(args$outcome) == 1) {
    if (args$outcome == "sens50") {
        knitr::kable(as.data.frame(weight_table), format = "latex", digits = 2,
                     caption = paste0("Table of Super Learner weights for each outcome, candidate learner ",
                                      "and cross-validation fold. We have removed `SL.' from the ",
                                      "name of each learner."),
                     label = "sl_weights", row.names = FALSE) %>%
            kableExtra::kable_styling(latex_options = c("scale_down")) %>%
            pack_rows("IC$_{50} < 1$", start_row = 11, end_row = 20, escape = FALSE) %>%
            cat(., file = here(ms_plots_dir, "supplement_sl_weight_table_sens50.tex"))
    } else if (args$outcome == "sens80") {
        knitr::kable(as.data.frame(weight_table), format = "latex", digits = 2,
                     caption = paste0("Table of Super Learner weights for each outcome, candidate learner ",
                                      "and cross-validation fold. We have removed `SL.' from the ",
                                      "name of each learner."),
                     label = "sl_weights", row.names = FALSE) %>%
            kableExtra::kable_styling(latex_options = c("scale_down")) %>%
            pack_rows("IC$_{80} < 1$", start_row = 1, end_row = 10, escape = FALSE) %>%
            cat(., file = here(ms_plots_dir, "supplement_sl_weight_table_sens80.tex"))
    } else {
        knitr::kable(as.data.frame(weight_table), format = "latex", digits = 2,
                     caption = paste0("Table of Super Learner weights for each outcome, candidate learner ",
                                      "and cross-validation fold. We have removed `SL.' from the ",
                                      "name of each learner."),
                     label = "sl_weights", row.names = FALSE) %>%
            kableExtra::kable_styling(latex_options = c("scale_down")) %>%
            pack_rows("IC$_{50}$ censored", start_row = 1, end_row = 10, escape = FALSE) %>%
            cat(., file = here(ms_plots_dir, "supplement_sl_weight_table_cens.tex"))
    }
} else {
    knitr::kable(as.data.frame(weight_table), format = "latex", digits = 2,
                 caption = paste0("Table of Super Learner weights for each outcome, candidate learner ",
                                  "and cross-validation fold. We have removed `SL.' from the ",
                                  "name of each learner."),
                 label = "sl_weights", row.names = FALSE) %>%
        kableExtra::kable_styling(latex_options = c("scale_down")) %>%
        pack_rows("IC$_{50}$ censored", start_row = 1, end_row = 10, escape = FALSE) %>%
        pack_rows("IC$_{50} < 1$", start_row = 11, end_row = 20, escape = FALSE) %>%
        cat(., file = here(ms_plots_dir, "supplement_sl_weight_table.tex"))
}

# plot of cv performance of the super learner in terms of AUC
cvauc_plots <- vector("list", length = length(args$outcome))
roc_plots <- vector("list", length = length(args$outcome))
nice_outcomes <- c(expression(bold(paste(IC[50], " Censored"))),
expression(bold(IC[50] < 1)), expression(bold(IC[80] < 1)))
all_cvaucs <- vector("list", length = length(args$outcome))
for (i in seq_len(length(args$outcome))) {
    o <- args$outcome[i]
    if (o == "cens") {
        o <- "ic50.censored"
    }
    fit_name <- paste0("cvfit_", o, ".rds")
    fit <- readRDS(here("code", "slfits", fit_name))
    class(fit) <- "myCV.SuperLearner"
    main_title <- nice_outcomes[i]
    cvauc_plots[[i]] <- plot(fit, main_title = main_title, xlim1 = 0.4, xlim2 = 1,
                             text_size = args$main_text_size, Rsquared = FALSE,
                             method = 'method.AUC')
    roc_plots[[i]] <- plot_roc_curves(fit, text_size = args$main_text_size,
                                      legend_text_size = args$axis_text_size * 0.8)
    all_cvaucs[[i]] <- cvAUC::ci.cvAUC(predictions = fit$SL.predict,
                                       labels = dat[!is.na(dat[, o]), o])
}
all_cvaucs
cvauc_plot <- plot_grid(plotlist = cvauc_plots, labels = "AUTO")
roc_plot <- plot_grid(plotlist = roc_plots, labels = "AUTO")
for (filetype in filetypes) {
    if (length(args$outcome) == 1) {
        perf_filename <- paste0("supplement_sl_perf_", args$outcome, ".")
        roc_filename <- paste0("supplement_sl_roc_", args$outcome, ".")
    } else {
        perf_filename <- paste0("supplement_sl_perf.")
        roc_filename <- paste0("supplement_sl_roc.")
    }
    ggsave(filename = here(ms_plots_dir, paste0(perf_filename, filetype)),
           plot = cvauc_plot, width = fig_width, height = fig_height,
           units = "in", dpi = 300)
    ggsave(filename = here(ms_plots_dir, paste0(roc_filename, filetype)),
           plot = roc_plot, width = fig_width, height = fig_height,
           units = "in", dpi = 300)
}
