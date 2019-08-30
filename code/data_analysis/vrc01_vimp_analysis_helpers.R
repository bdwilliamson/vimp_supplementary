## helper functions for VRC01 data analysis
# match the correct outcome to a given row in the table
match_y <- function(row, y1, y2, folds1, folds2, ord =  c("ic50", "ic80", "slope_mod", "sens.resis", "cens")) {
  if (as.numeric(row[3]) == 1) {
    tmp_y <- y1[[which(as.character(row[1]) == ord)]]  
    # ret_y <- make_y_lst(tmp_y, folds1[[which(as.character(row[1]) == ord)]])
    ret_y <- list(y = tmp_y, folds = folds1[[which(as.character(row[1]) == ord)]])
  } else {
    tmp_y <- y2[[which(as.character(row[1]) == ord)]]
    # ret_y <- make_y_lst(tmp_y, folds2[[which(as.character(row[1]) == ord)]])
    ret_y <- list(y = tmp_y, folds = folds2[[which(as.character(row[1]) == ord)]])
  }
  return(ret_y)
}

# match an outcome, e.g. "Y", to a character value, e.g. "ic50"
match_chr_y <- function(o, s, outs = c("Y", "Y2", "Y.80", "Y2.80", "Y.slope", "Y2.slope", "Y.sens.resis", "Y2.sens.resis", "Y.cens", "Y2.cens"),
 ord = c("ic50", "ic80", "slope_mod", "sens.resis", "cens")) {
    # get which one it is
    indx <- which(o == ord)
    # make a vector of repeats to choose from based on o and s
    reps <- rep(1:length(unique(grp_mat$o)), each = 2)
    # return the one based on o and s
    ret_indx <- which(indx == reps)[s]
    ret <- outs[ret_indx]
    return(ret)
}

# turn each vector in the y list into a list, based on V folds
make_y_lst <- function(y, folds) {
  ## break it up into a list
  y_lst <- list()
  flds <- sort(unique(folds))
  for (v in 1:max(flds)) {
    y_lst[[v]] <- y[folds == flds[v]]
  }
  return(y_lst)
}
# create a giant matrix with results
# and alphabetize the groups, then the individual features

create_table <- function(nms_grp, nms_ind, p_nms_grp, p_nms_ind, g_mat, i_mat, col_nms, num_measures = 5, num_datasets = 2, alphabetize = TRUE) {
  # first make a matrix with the group ones
  grp <- matrix(0, nrow = length(p_nms_grp), ncol = num_measures*length(unique(grp_mat$o))*num_datasets)
  # for each feature set, create the matrix
  for (i in 1:length(p_nms_grp)) {
    # get the feature set names
    f_set_nms <- get_feature_list(nms_grp, i, g_mat)
    # pass them in to get the row
    grp[i, ] <- make_matrix(f_set_nms, g_mat, num_measures)
  }
  # now slap on row names
  rownames(grp) <- p_nms_grp
  # alphabetize
  if (alphabetize) {
    ret_grp <- grp[order(rownames(grp)), ]  
  } else { # do it by order of the estimate for cens
    ret_grp <- grp[order(grp[, 1], decreasing = TRUE), ]
  }
  
  
  # do the same thing for independent
  ind <- matrix(0, nrow = length(p_nms_ind), ncol = num_measures*length(unique(grp_mat$o))*num_datasets)
  for (i in 1:length(p_nms_ind)) {
    # get feature set names
    f_set_nms <- get_feature_list(nms_ind, i, i_mat)
    # make a row
    ind[i, ] <- make_matrix(f_set_nms, i_mat, num_measures)
  }
  # row names
  rownames(ind) <- p_nms_ind
  # alphabetize
  if (alphabetize) {
    ret_ind <- ind[order(rownames(ind)), ]  
  } else {
    ret_ind <- ind[order(ind[, 1], decreasing = TRUE), ]  
  }
  
  
  # make a big matrix
  ret <- rbind(ret_grp, ret_ind)
  colnames(ret) <- col_nms
  return(ret)
}

# helper function that takes in a name and returns the estimate, SE, and CI formatted nicely
return_est_ci <- function(nm, num_measures = 4) {
  if (num_measures == 4) {
    eval(parse(text = paste0("ret <- c(", nm, "$est, ", nm, "$se, ", nm, "$ci)")))
  } else if (num_measures == 5) {
    eval(parse(text = paste0("ret <- c(", nm, "$est, ", nm, "$se, ", nm, "$ci, ", nm,  "$test)")))
  } else {
    eval(parse(text = paste0("ret <- c(", nm, "$est, ", nm, "$se, ", nm, "$ci, ", nm,  "$test, ", nm, "$p_value)")))
  }
  
  return(ret)
}

# helper function that takes in a list of names and makes a matrix with columns
# est, ci for each outcome
make_matrix <- function(lst, mat, num_measures) {
  # get the estimates for each name in the list, combine it into a vector
  tmp <- lapply(lst, return_est_ci, num_measures)
  ret <- do.call(c, tmp)
  return(ret)
}

# helper function to extract a list the desired feature/feature set from a list of names
get_feature_list <- function(lst, num, mat) {
  indx <- which(mat$f == num)
  ret <-  as.list(lst[indx])
}
# create a table of results for each outcome, with the average as well
create_sub_table <- function(avg, set1, set2, nms, num_measures) {
  # order set 1 and set 2 based on avg
  ord_1 <- match(avg$s, set1$s)
  ord_2 <- match(avg$s, set2$s)
  
  ## if any are NA, set to the end
  if (any(is.na(ord_1))) ord_1[is.na(ord_1)] <- length(avg$s)
  if (any(is.na(ord_2))) ord_2[is.na(ord_2)] <- length(avg$s)
  
  # now combine them together based on the matching
  ret <- cbind(avg$mat[, 1:num_measures], set1$mat[ord_1, 1:num_measures], set2$mat[ord_2, 1:num_measures])
  if (num_measures == 4) {
    colnames(ret) <- c("avg_est", "avg_se", "avg_cil", "avg_ciu",
                       "1_est", "1_se", "1_cil", "1_ciu",
                       "2_est", "2_se", "2_cil", "2_ciu")  
  } else if (num_measures == 5) {
    colnames(ret) <- c("avg_est", "avg_se", "avg_cil", "avg_ciu", "avg_test",
                       "1_est", "1_se", "1_cil", "1_ciu", "1_test",
                       "2_est", "2_se", "2_cil", "2_ciu", "2_test")
  } else {
    colnames(ret) <- c("avg_est", "avg_se", "avg_cil", "avg_ciu", "avg_test", "avg_p_value",
                       "1_est", "1_se", "1_cil", "1_ciu", "1_test", "1_p_value",
                       "2_est", "2_se", "2_cil", "2_ciu", "2_test", "2_p_value")
  }
  
  rownames(ret) <- nms[as.numeric(avg$s)]
  return(ret)
}

create_outcome_table <- function(avg_grp, set1_grp, set2_grp, nms_grp, avg_ind, set1_ind, set2_ind, nms_ind, col_nms, num_measures) {
  # sub table for groups
  grp <- create_sub_table(avg_grp, set1_grp, set2_grp, nms_grp, num_measures)
  # sub table for individual
  ind <- create_sub_table(avg_ind, set1_ind, set2_ind, nms_ind, num_measures)
  # row bind together
  ret <- rbind(grp, ind)
  colnames(ret) <- col_nms
  return(ret)
}
