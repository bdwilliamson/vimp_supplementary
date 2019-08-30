##########################################################################
##
## FILE: vrc01_vimp_analysis.R
##
## PURPOSE: variable importance analysis of binary ic50 endpoints
##          from VRC01 CATNAP analysis, using deviance, accuracy, auc
##
##########################################################################

##------------------------------------------------------------------------
## set up
##------------------------------------------------------------------------
## load required functions and packages
# only run this if necessary to update package
# install.packages("~/Projects/UW/vimp_1.3.0.tar.gz", repos = NULL, type = "source")
library("vimp")
library("argparse")
code_dir <- "../code/"
plots_dir <- "../code/plots/"
results_dir <- "../code/vrc01_analysis_output/"
source(paste0(code_dir, "vrc01_vimp_analysis_helpers.R"))

parser <- ArgumentParser()
parser$add_argument("--vimp-measure", default = "deviance", help = "the variable importance measure of interest")
parser$add_argument("--level", type = "double", default = 0.95, help = "nominal level for confidence intervals")
parser$add_argument("--point-size", type = "double", default = 2, help = "point size")
parser$add_argument("--axis-text-size", type = "double", default = 2.5, help = "axis text size")
parser$add_argument("--plot-text-size", type = "double", default = 2, help = "plot text size")
parser$add_argument("--main-text-size", type = "double", default = 3, help = "main text size")
args <- parser$parse_args()

## set directories, dataset version, etc.
vrc01_results_dir <- "~/Projects/VIDD/hvtn/vrc01_genotypic_resistance_score/results/"
vrc01_code_dir <- "~/Projects/VIDD/hvtn/vrc01_genotypic_resistance_score/code/"
vrsn <- "12"
nms_pub <- list(c(expression("Sensitive/Resistant Only; deviance"),
             expression(paste(IC[50], " Censored; deviance", sep = ""))),
             c(expression("Sensitive/Resistant Only; accuracy"),
             expression(paste(IC[50], " Censored; accuracy", sep = ""))),
             c(expression("Sensitive/Resistant Only; AUC"),
             expression(paste(IC[50], " Censored; AUC", sep = ""))))
run <- "standard"

nms_group <- c("VRC01 contact sites (Set 1)", "CD4 binding sites (Set 2)", "ESA (Set 3)",
                "Glycosylation sites (Set 4)", "Covarying sites (Set 5)", "PNG sites (Set 6)",
                "gp41 sites (Set 7)",
                "N-linked glycosylation (Set 8)", "Subtype (Set 9)", "Sequons (Set 10)",
                "Viral geometry (Set 11)", "Cysteines (Set 12)", "Steric bulk (Set 13)")
num_indiv <- 797
num_grp <- 13

## set up the groups
grp_mat <- expand.grid(o = c("sens.resis", "cens"),
                       f = 1:num_grp, s = 1:2, stringsAsFactors = FALSE)
ind_mat <- expand.grid(o = c("sens.resis", "cens"),
                       f = 1:num_indiv, s = 1:2, stringsAsFactors = FALSE)

## -----------------------------------------------------------------------
## variable importance for each outcome
## -----------------------------------------------------------------------
## only run this if I need to re-compile results
# source(paste0(code_dir, "vrc01_vimp_analysis_compute_group_vimp.R"))
# source(paste0(code_dir, "vrc01_vimp_analysis_compute_indiv_vimp.R"))
## read in group results to pass to plot
for (i in 1:length(unique(grp_mat$o))) {
    eval(parse(text = paste0("grp_", unique(grp_mat$o)[i], "_avg <- readRDS(file = paste0(results_dir, args$vimp_measure, '_vim_avg_grp_", unique(grp_mat$o)[i], ".rds'))")))
    for (j in 1:2) {
        eval(parse(text = paste0("grp_", unique(grp_mat$o)[i], "_", j, " <- readRDS(file = paste0(results_dir, args$vimp_measure, '_vim_grp_", unique(grp_mat$o)[i], "_", j, ".rds'))")))
    }
}


##------------------------------------------------------------------------
## make nice plots
##------------------------------------------------------------------------
fig_width <- 2590
fig_height <- fig_width

## plot that combines results for same outcome into one figure (makes 5 total figures)
## each plot has results for dataset 1 (red circles) and dataset 2 (blue diamonds),
## along with the average across the two (black squares), with CIs
for (i in 1:length(unique(grp_mat$o))) {
  # names
  eval(parse(text = paste0("nms_tmp <- unique(grp_mat$f)[as.numeric(unlist(lapply(strsplit(rownames(grp_", unique(grp_mat$o)[i], "_avg$mat), '_', fixed = TRUE), function(x) x[4])))]")))
  nms <- rev(nms_tmp)
  corrections <- c(-0.2, -0.4)
  # plot it
  tmp1 <- eval(parse(text = paste0("grp_", unique(grp_mat$o)[i], "_avg")))
  ## add text with CI values
  tmp_cis <- round(tmp1$mat[, c("cil", "ciu")], 3)
  text_cis <- apply(tmp_cis, 1, function(x) paste0("[", x[1], ", ", x[2], "]"))

  # make a png
  png(paste0(plots_dir, "vim_grp_combined_outcome_", unique(grp_mat$o)[i], "_", run, "_level_", args$level, "_", args$vimp_measure, ".png"), width = fig_width, height = fig_height, units = "px", res = 300)
  plot(tmp1, nms, mar = c(5.1, 6, 4, 1), col = 'black', pch = 15, cex.axis = args$axis_text_size, cex = args$point_size, xlim = c(0, 0.4),
       yaxt = "n", xaxt = "n")
  axis(side = 2, at = 1:length(tmp1$mat$est) + corrections[1], labels = nms, cex.axis = args$axis_text_size, tick = FALSE, las = 2)
  axis(side = 1, at = seq(0, 0.4, 0.05), line = 0, cex.axis = args$axis_text_size, mgp = c(3, 1.5, 1))
  ## add left-hand arrowheads for ones that go outside plot margins
  outside_left_margin <- tmp1$mat$cil < -0.015
  outside_right_margin <- tmp1$mat$ciu > 0.415
  points(rep(-0.015, length(tmp1$mat$cil[outside_left_margin])), rev(1:length(tmp1$mat$cil))[outside_left_margin], pch = -9668)
  points(rep(0.415, length(tmp1$mat$cil[outside_right_margin])), rev(1:length(tmp1$mat$cil))[outside_right_margin], pch = -9658)
  text(rep(0.35, length(tmp1$mat$ciu)), rev(1:length(tmp1$mat$ciu)), text_cis, cex = args$plot_text_size)
  # plot the individual ones
  cols <- c("red", "blue")
  pchs <- c(16, 18)
  for (j in 1:2) {
    tmp <- eval(parse(text = paste0("grp_", unique(grp_mat$o)[i], "_", j)))
    ord_tmp <- tmp$mat[match(as.numeric(tmp1$s), tmp$s), ]
    ord_tmp_2 <- ord_tmp[nrow(ord_tmp):1, ]
    points(ord_tmp_2$est, 1:num_grp + corrections[j], col = cols[j], pch = pchs[j], cex = args$point_size)
    arrows(unlist(ord_tmp_2$cil), 1:dim(ord_tmp_2)[1] + corrections[j], unlist(ord_tmp_2$ciu), 1:dim(ord_tmp_2)[1] + corrections[j],
           length = 0, angle = 90, lwd = 2, col = cols[j])
    outside_left_margin <- tmp$mat$cil < -0.015
    outside_right_margin <- tmp$mat$ciu > 0.415
    points(rep(-0.015, length(tmp1$mat$cil[outside_left_margin])), rev(1:length(tmp$mat$cil))[outside_left_margin] + corrections[j], pch = -9668, col = cols[j])
    points(rep(0.415, length(tmp1$mat$cil[outside_right_margin])), rev(1:length(tmp$mat$cil))[outside_right_margin] + corrections[j], pch = -9658, col = cols[j])
  }
  # add labels
  title(main = unlist(nms_pub[args$vimp_measure == c("deviance", "accuracy", "auc")])[i], cex.main = args$main_text_size, font.main = 1)
  # bump out the y-axis label, slightly
  # title(ylab = "Feature Group", cex.lab = args$axis_text_size)
  mtext("Feature Group", side = 2, line = 4, cex = args$axis_text_size)
  title(xlab = "Estimated Variable Importance", cex.lab = args$axis_text_size, line = 4)
  ## add a legend, only for the first one
  if (args$vimp_measure == "accuracy") {
    legend("bottom", legend = c("Average", "Dataset 1", "Dataset 2"), col = c("black", "red", "blue"),
           pch = c(15, pchs), cex = args$axis_text_size)
  }
  # close it off
  dev.off()

  # make a tiff
  tiff(paste0(plots_dir, "vim_grp_combined_outcome_", unique(grp_mat$o)[i], "_", run, "_level_", args$level, "_", args$vimp_measure, ".tiff"), width = fig_width, height = fig_height, units = "px", res = 300)
  plot(tmp1, nms, mar = c(5.1, 6, 4, 1), col = 'black', pch = 15, cex.axis = args$axis_text_size, cex = args$point_size, xlim = c(0, 0.4),
       yaxt = "n", xaxt = "n")
  axis(side = 2, at = 1:length(tmp1$mat$est) + corrections[1], labels = nms, cex.axis = args$axis_text_size, tick = FALSE, las = 2)
  axis(side = 1, at = seq(0, 0.4, 0.05), line = 0, cex.axis = args$axis_text_size, mgp = c(3, 1.5, 1))
  ## add left-hand arrowheads for ones that go outside plot margins
  outside_left_margin <- tmp1$mat$cil < -0.015
  outside_right_margin <- tmp1$mat$ciu > 0.415
  points(rep(-0.015, length(tmp1$mat$cil[outside_left_margin])), rev(1:length(tmp1$mat$cil))[outside_left_margin], pch = -9668)
  points(rep(0.415, length(tmp1$mat$cil[outside_right_margin])), rev(1:length(tmp1$mat$cil))[outside_right_margin], pch = -9658)
  text(rep(0.35, length(tmp1$mat$ciu)), rev(1:length(tmp1$mat$ciu)), text_cis, cex = args$plot_text_size)
  # plot the individual ones
  cols <- c("red", "blue")
  pchs <- c(16, 18)
  for (j in 1:2) {
    tmp <- eval(parse(text = paste0("grp_", unique(grp_mat$o)[i], "_", j)))
    ord_tmp <- tmp$mat[match(as.numeric(tmp1$s), tmp$s), ]
    ord_tmp_2 <- ord_tmp[nrow(ord_tmp):1, ]
    points(ord_tmp_2$est, 1:num_grp + corrections[j], col = cols[j], pch = pchs[j], cex = args$point_size)
    arrows(unlist(ord_tmp_2$cil), 1:dim(ord_tmp_2)[1] + corrections[j], unlist(ord_tmp_2$ciu), 1:dim(ord_tmp_2)[1] + corrections[j],
           length = 0, angle = 90, lwd = 2, col = cols[j])
    outside_left_margin <- tmp$mat$cil < -0.015
    outside_right_margin <- tmp$mat$ciu > 0.415
    points(rep(-0.015, length(tmp1$mat$cil[outside_left_margin])), rev(1:length(tmp$mat$cil))[outside_left_margin] + corrections[j], pch = -9668, col = cols[j])
    points(rep(0.415, length(tmp1$mat$cil[outside_right_margin])), rev(1:length(tmp$mat$cil))[outside_right_margin] + corrections[j], pch = -9658, col = cols[j])
  }
  # add labels
  title(main = unlist(nms_pub[args$vimp_measure == c("deviance", "accuracy", "auc")])[i], cex.main = args$main_text_size, font.main = 1)
  # bump out the y-axis label, slightly
  # title(ylab = "Feature Group", cex.lab = args$axis_text_size)
  mtext("Feature Group", side = 2, line = 4, cex = args$axis_text_size)
  title(xlab = "Estimated Variable Importance", cex.lab = args$axis_text_size, line = 4)
  ## add a legend, only for the first one
  if (args$vimp_measure == "accuracy") {
    legend("bottom", legend = c("Average", "Dataset 1", "Dataset 2"), col = c("black", "red", "blue"),
           pch = c(15, pchs), cex = args$axis_text_size)
  }
  ## close it off
  dev.off()
}

##------------------------------------------------------------------------
## make nice tables
##------------------------------------------------------------------------
# make the column names
measures <- c("Est.", "SE", "CIL", "CIU", "VIMP > 0", "p_value")
num_measures <- length(measures)
col_nms <- paste0("Dataset ", rep(1:2, each = num_measures*length(unique(full_mat$o))), ": ", rep(rep(c("cens", "sens.resis"), each = num_measures), 2), ";", rep(measures, length(unique(full_mat$o))*2))

# make a giant matrix with estimate, ci for each outcome in columns and feature (set) in rows
mat <- create_table(grp_nms, ind_nms, nms_group, nms_ind, grp_mat, ind_mat, col_nms, num_measures = num_measures, num_datasets = 2, TRUE)

# save off the giant table
write.csv(mat, paste0(plots_dir, "vimp_table_with_groups_", run, "_level_", args$level, "_", args$vimp_measure, ".csv"))

# now make a giant matrix ordered by importance rather than by alphabetical
mat2 <- create_table(grp_nms, ind_nms, nms_group, nms_ind, grp_mat, ind_mat, col_nms, num_measures = num_measures, num_datasets = 2, FALSE)

# save off the giant table
write.csv(mat2, paste0(plots_dir, "vimp_table_with_groups_by_importance_", run, "_level_", args$level, "_", args$vimp_measure, ".csv"))

# now make a table for each outcome
# ordered by the average (groups first, then individual)
# average, set 1, set 2
# can choose a cutoff for publication ready table
col_nms_2 <- paste0(rep(c("Avg: ", "Data set 1: ", "Data set 2: "), each = num_measures), rep(measures, 3))
for (i in 1:length(unique(full_mat$o))) {
    eval(parse(text = paste0(unique(full_mat$o)[i], " <- create_outcome_table(grp_", unique(full_mat$o)[i], "_avg",
                             ", grp_", unique(full_mat$o)[i], "_1", ", grp_", unique(full_mat$o)[i], "_2",
                             ", nms_group , ind_", unique(full_mat$o)[i], "_avg", ", ind_", unique(full_mat$o)[i], "_1",
                             ", ind_", unique(full_mat$o)[i], "_2", ", nms_ind, col_nms_2, num_measures)")))

  eval(parse(text = paste0("write.csv(", unique(full_mat$o)[i], ", paste0(plots_dir, 'vimp_table_by_outcome_", unique(full_mat$o)[i], "_", run, "_level_", args$level, "_", args$vimp_measure, ".csv'))")))
}
