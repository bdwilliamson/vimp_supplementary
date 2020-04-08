#!/usr/local/bin/Rscript

## define variable groups for Magaret et al. (2019) paper

## function to set a variable group
get_variable_group <- function(pred_names, hxb2_sites) {
    ## get only the sites from the pred_names vector
    aa_positions <- unlist(lapply(strsplit(pred_names, ".", fixed = TRUE), function(x) x[2]))

    ## get the ones matching the input sites
    site_character_vars <- pred_names[aa_positions %in% hxb2_sites]
    return(site_character_vars)
}

get_variable_groups <- function(data) {
    ## set up variable groups
    ## geographic region variables
    geog_vars <- colnames(dat)[grepl("geographic.region.of.origin.is", colnames(dat))]
    ## subtype variables
    subtype_vars <- colnames(dat)[grepl("subtype.is", colnames(dat))]
    ## Group 1: VRC01 contact sites
    vrc01_contact_sites <- c(97, 123, 124, 198, 276, 278, 279, 280, 281, 282, 365, 366, 367, 368, 371, 427, 428, 429, 430, 455, 456, 457, 458, 459, 460, 461, 463, 465, 466, 467, 469, 472, 473, 474, 476)
    ## Group 2: CD4 binding sites
    cd4_binding_sites <- c(124, 125, 126, 127, 196, 198, 279, 280, 281, 282, 283, 365, 366, 367, 368, 369, 370, 374, 425, 426, 427, 428, 429, 430, 431, 432, 455, 456, 457, 458, 459, 460, 461, 469, 471, 472, 473, 474, 475, 476, 477)
    ## Group 3: ESA sites
    esa_sites <- c(97, 198, 276, 278, 279, 280, 281, 282, 365, 366, 367, 368, 371, 415, 428, 429, 430, 455, 457, 458, 459, 460, 461, 467, 469, 473, 474, 476)
    ## Group 4: glyco sites
    glyco_sites <- c(61, 64, 197, 276, 362, 363, 386, 392, 462, 463)
    ## Group 5: covarying sites
    covar_sites <- c(46, 132, 138, 144, 150, 179, 181, 186, 190, 290, 321, 328, 354, 389, 394, 396, 397, 406)
    ## Group 6: vrc01-specific PNG sites
    png_sites <- c(130, 139, 143, 156, 187, 197, 241, 262, 289, 339, 355, 363, 406, 408, 410, 442, 448, 460, 462)
    ## Group 7: gp41 sites
    gp41_sites <- c(544, 569, 582, 589, 655, 668, 675, 677, 680, 681, 683, 688, 702)
    ## Group 8: potential n-linked glycosylation sites not in VRC01 contact sites or paratope or sites with covariability
    pngs_ind_1 <- which(colnames(data) == "hxb2.29.sequon_actual.1mer")  # the first index
    pngs_ind_2 <- which(colnames(data) == "hxb2.824.sequon_actual.1mer")  # the last index
    pngs_init <- colnames(data)[pngs_ind_1:pngs_ind_2]
    pngs_aa_pos <- unlist(lapply(strsplit(pngs_init, ".", fixed = TRUE), function(x) x[2]))
    pngs_novrc01_sites <- pngs_aa_pos[!(pngs_aa_pos %in% vrc01_contact_sites) & !(pngs_aa_pos %in% sort(c(cd4_binding_sites, esa_sites, glyco_sites, covar_sites, png_sites)))]
    ## Group 9: majority virus subtypes (already created, subtype_vars above)
    ## Group 10: region-specific counts of PNG sites
    glycosylation_vars <- colnames(dat)[grepl("sequons", colnames(dat))]
    ## Group 11: viral geometry
    geometry_vars <- colnames(dat)[grepl("length", colnames(dat))]
    ## Group 12: cysteine counts
    cysteine_vars <- colnames(dat)[grepl("cysteines", colnames(dat))]
    ## Group 13: steric bulk at critical locations
    steric_bulk_vars <- colnames(dat)[grepl("taylor", colnames(dat))]

    ## get all variable groups
    all_data_names <- colnames(data)
    aa_vrc01_vars <- get_variable_group(all_data_names, vrc01_contact_sites)
    aa_cd4bs_vars <- get_variable_group(all_data_names, cd4_binding_sites)
    aa_esa_vars <- get_variable_group(all_data_names, esa_sites)
    aa_glyco_vars <- get_variable_group(all_data_names, glyco_sites)
    aa_covar_vars <- get_variable_group(all_data_names, covar_sites)
    aa_png_vars <- get_variable_group(all_data_names, png_sites)
    aa_gp41_vars <- get_variable_group(all_data_names, gp41_sites)
    aa_pngs_novrc01_vars <- get_variable_group(all_data_names, pngs_novrc01_sites)
    return(list(vrc01 = aa_vrc01_vars, cd4bs = aa_cd4bs_vars, esa = aa_esa_vars, glyco = aa_glyco_vars, covar = aa_covar_vars, pngs = aa_png_vars, gp41 = aa_gp41_vars, pngs_novrc01 = aa_pngs_novrc01_vars, subtype = subtype_vars, sequons = glycosylation_vars, geometry = geometry_vars, cysteines = cysteine_vars, steric_bulk = steric_bulk_vars, geog = geog_vars))
}
