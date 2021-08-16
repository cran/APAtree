#' Calculate the index for neighborhood diversity NDiv for an apa_list
#'
#' NDiv assesses the average dissimilarity between trees and their neighbors
#' using APA-maps. Dissimilarity can be assessed using a species- or a
#' trait-based approach. Upscaling of the tree-level NDiv to the species- or
#' stand-level are provided.
#'
#' @template apa_list_arg
#'
#' @inheritParams apa_list
#'
#' @param pdiv If TRUE (default), pdiv will be calculated (proportion-based diversity,
#'   average dissimilarity between a tree and all other trees in a stand,
#'   irrespective of the spatial configuration.)
#'   
#' @param dis_trait_column A list containing combinations of traits that will
#'   be used to estimate dissimilarity between trees when calculating
#'   apa-properties. Refers to column names of \code{tree_dat}.
#' 
#' @details See Glatthorn (2021) for details.
#' 
#' @template ref_glatthorn_2021
#' 
#' @example examples/example_apa_ndiv.R
#' 
#' @export
#' 

apa_ndiv <- function(apa_list, 
                     dis_trait_column, 
                     dis_method = "gowdis",
                     dis_transform = sqrt, 
                     scope = "global",
                     edge_correction = "none",
                     pdiv = TRUE){
  apa_config <- attr(apa_list, "apa_config")
  if(is.null(dis_trait_column)){
    stop("A `dis_trait_column` has to be specified for the calculation of `ndiv`.")
  }
  if(is.null(edge_correction)){
    edge_correction <- apa_config$edge_correction
  }
  if(!is.na(apa_config$edge_correction) && apa_config$edge_correction != edge_correction){
    stop("There can only be a single edge correction method for one 'apa_list' object.")
  }
  if(!"apa_size" %in% apa_config$apa_properties){
    apa_list <- apa_size(apa_list, edge_correction = edge_correction)
    apa_config <- attr(apa_list, "apa_config")
  }
  if(any(grepl("ndiv", apa_config$apa_properties))){
    stop("`apa_ndiv` has already been calculated. Please use `apa_drop_properties()` first to recalculate.")
  }
  if(! scope %in% c("global", "local")){
    stop("`scope` has to be either `global` or `local`.")
  }
  apa_map_tree_layer <- 
    lapply(apa_list$plot_dat$apa_map, raster::subset, apa_config$tree_id_column)
  apa_map_tree_levels <- 
    lapply(apa_map_tree_layer, raster::levels)
  apa_map_tree_count <- c(lapply_deep(apa_map_tree_levels, nrow), recursive = TRUE)
  if(any(apa_map_tree_count == 1)){
    #stop("There is only a single tree in apa-map `", paste(names(apa_map_tree_count), collapse = "`, `"),
    #     "`. Calculation of ndiv requires at least two trees.\n")
  }
  if(!is.list(dis_trait_column)){
    dis_trait_column <- list(dis_trait_column)
  }
  if(is.null(names(dis_trait_column))){
    dist_trait_lengths <- sapply(dis_trait_column, length)
    multi_trait_names <- dis_trait_column[dist_trait_lengths > 1]
    multi_trait_names <- lapply(multi_trait_names, substr, 1, 1)
    multi_trait_names <- lapply(multi_trait_names, paste, collapse = "")
    trait_names <- dis_trait_column
    trait_names[dist_trait_lengths > 1] <- multi_trait_names
    names(dis_trait_column) <- trait_names
  }
  if(any(duplicated(names(dis_trait_column)))){
    stop("Please choose unique names for `dis_trait_column`.")
  }
  dis_trait_missing <- 
    lapply(dis_trait_column, setdiff, y = names(apa_list$tree_dat))
  dis_trait_missing <- c(dis_trait_missing, recursive = TRUE)
  if(length(dis_trait_missing) > 0){
    stop("Columns `", paste(dis_trait_missing, collapse = "`, `"), "` are missing in `tree_dat`.")
  }
  apa_config$dis_trait_column <- dis_trait_column
  add_to_apa_properties <- paste0(names(dis_trait_column), "_ndiv")
  if(pdiv){
    add_to_apa_properties <- 
      c(add_to_apa_properties, 
        paste0(names(dis_trait_column), "_pdiv"))
  }
  apa_config$apa_properties <- 
    as.vector(stats::na.omit(unique(c(apa_config$apa_properties, add_to_apa_properties))))
  apa_config$dis_method <- dis_method
  apa_config$dis_transform <- dis_transform
  apa_config$scope <- scope
  apa_map <- apa_list$plot_dat$apa_map
  trait_dat <- sf::st_drop_geometry(apa_list$tree_dat)
  trait_dat <- lapply(lapply(apa_config$dis_trait_column, c, apa_config$tree_id_column),
                      drop_columns_except,
                      x = trait_dat)
  if(apa_config$core_column == apa_config$buffer_column){
    top_level_plot <- NULL
    top_level_id_column <- NULL
  }else{
    top_level_plot <- apa_list$plot_dat
    sf::st_geometry(top_level_plot) <- apa_config$core_column
    top_level_plot <- top_level_plot[, apa_config$plot_id_column, drop = FALSE]
    top_level_id_column <- c(plot_dat = apa_config$plot_id_column)
  }
  border_tree_id <- 
    stats::na.omit(apa_list$tree_dat[[apa_config$tree_id_column]][apa_list$tree_dat$border_tree])
  output <- NULL
  for(trait_i in names(trait_dat)){
    message(paste0("\nAggregating `", trait_i, "` neighborhood diversity:"), appendLF = FALSE)
    output_i <- apa_ndiv_calc(unclass(apa_list), apa_config, apa_map,
                              trait_dat = trait_dat[[trait_i]], subplot = top_level_plot,
                              subplot_id_column = top_level_id_column, pdiv = pdiv,
                              edge_correction = edge_correction, border_tree_id = border_tree_id)
    subplot_list <- apa_list$subplot_dat
    for(subplot_i in names(apa_config$subplot_id_column)){
      message(paste0("`", subplot_i, "` - aggregating neighborhood diversity:  "), appendLF = FALSE)
      subplot_dat_i <- apa_list$subplot_dat[[subplot_i]]
      subplot_list[[subplot_i]] <- 
        apa_ndiv_calc(subplot_dat_i, apa_config, apa_map,
                      trait_dat = trait_dat[[trait_i]],
                      subplot = subplot_dat_i[[subplot_i]],
                      subplot_id_column = apa_config$subplot_id_column[subplot_i],
                      pdiv = pdiv, edge_correction = edge_correction,
                      border_tree_id = border_tree_id)
    }
    message("")
    output_i$subplot_dat <- subplot_list
    output_i <- 
      lapply_deep(output_i, .what = "data.frame",
                  function(x){
                    names(x)[names(x) == "ndiv"] <- paste0(trait_i, "_", "ndiv")
                    names(x)[names(x) == "pdiv"] <- paste0(trait_i, "_", "pdiv")
                    x
                  })
    if(is.null(output)){
      output <- output_i
    }else{
      output <- mapply_deep(list(output, output_i), .what = "data.frame",
                            .f = add_col, col = paste0(trait_i, "_ndiv"))
      output <- mapply_deep(list(output, output_i), .what = "data.frame",
                            .f = add_col, col = paste0(trait_i, "_pdiv"))
    }
  }
  output <- new_apa_list(output, apa_config = apa_config)
  output
}

#' @keywords internal
apa_ndiv_calc <- 
  function(dat_list, apa_config, apa_map, trait_dat, subplot = NULL,
           subplot_id_column = NULL, pdiv = TRUE, prefix = NULL,
           edge_correction = "none", border_tree_id = NULL){
    if(!is.na(apa_config$edge_correction) && apa_config$edge_correction != edge_correction){
      stop("There can only be a single edge correction method for one 'dat_list' object.")
    }
    if(!edge_correction[1] %in% c("none", "critical", "border_tree_exclusion")){
      stop("`edge_correction` has to be one of `none`, `critical` or `border_tree_exclusion`.")
    }
    if(edge_correction == "critical"){
      critical_layer <- "critical"
    }else{
      critical_layer <- NULL
    }
    
    
    # Check if observations in boundary_length_dat are present in trait_dat
    if(inherits(trait_dat, "sf")){
      trait_dat <- sf::st_drop_geometry(trait_dat)
    }
    if(is.null(subplot)){
      agg_column <- apa_config$plot_id_column
      
    }else{
      agg_column <- subplot_id_column
    }
    if(edge_correction != "border_tree_exclusion" | is.null(subplot)){
      boundary_length_dat <- 
        boundary_length(rst = apa_map,
                        subplot = subplot,
                        plot_id_column = apa_config$plot_id_column,
                        critical_layer = critical_layer,
                        remove_na = FALSE)
    }else{
      boundary_length_dat_values <- 
        boundary_length(rst = apa_map,
                        subplot = NULL, # border tree correction uses boundary length outside the core zone 
                        plot_id_column = apa_config$plot_id_column,
                        critical_layer = critical_layer,
                        remove_na = FALSE)
      boundary_length_dat <- sf::st_join(subplot, dat_list$tree_dat[, apa_config$tree_id_column])
      boundary_length_dat <- sf::st_drop_geometry(boundary_length_dat)
      sfc_idx <- sapply(boundary_length_dat, inherits, "sfc")
      boundary_length_dat <- boundary_length_dat[!sfc_idx]
      boundary_length_dat <- merge(boundary_length_dat, boundary_length_dat_values,
                                   by = c(apa_config$plot_id_column, apa_config$tree_id_column))
    }
    if(edge_correction == "critical"){
      boundary_length_dat$boundary_length <- boundary_length_dat$boundary_length * (1- boundary_length_dat[[critical_layer]])
    } else if(edge_correction == "border_tree_exclusion"){
      border_tree_idx <- boundary_length_dat[[apa_config$tree_id_column]] %in% border_tree_id
      border_tree_neighbor_idx <- 
        boundary_length_dat[[paste0(apa_config$tree_id_column, "_bc")]] %in% border_tree_id
      boundary_length_dat$border_tree <- border_tree_idx | border_tree_neighbor_idx
      if(any(is.na(boundary_length_dat[[paste0(apa_config$tree_id_column, "_bc")]]) & !boundary_length_dat$border_tree)){
        stop("BUUUH!")
      }
      boundary_length_dat$boundary_length <- boundary_length_dat$boundary_length * (1- boundary_length_dat$border_tree)
    }
    
    trait_dat_trees <- trait_dat[[apa_config$tree_id_column]]
    boundary_length_dat_trees <- boundary_length_dat[, paste0(apa_config$tree_id_column, c("", "_bc"))]
    boundary_length_dat_trees <- unique(c(boundary_length_dat_trees, recursive = TRUE))
    missing_trees <- stats::na.omit(setdiff(boundary_length_dat_trees, trait_dat_trees))
    if(length(missing_trees) > 0){
      stop("Observations '", paste(missing_trees, collapse = "', '"), "' are missing in 'trait_dat'.")
    }
    tree_id <- trait_dat[[apa_config$tree_id_column]]
    if(length(unique(tree_id)) != length(tree_id) | any(is.na(tree_id))){
      stop("'tree_id_column' has to point to a unique identifier without 'NA' values.")
    }
    trait_val <- as.data.frame(trait_dat)[setdiff(names(trait_dat), apa_config$tree_id_column)]
    fac_trait_idx <- sapply(trait_val, class) == "factor"
    trait_val[fac_trait_idx] <- lapply(trait_val[fac_trait_idx], as.character)
    rownames(trait_val) <- trait_dat[[apa_config$tree_id_column]]
    if(apa_config$scope == "global"){
      minmax_dummies <- 
        rbind(as.data.frame(lapply(trait_val, min)),
              as.data.frame(lapply(trait_val, max)))
      rownames(minmax_dummies) <- c(".min_dummy", ".max_dummy")
    }
    
    if(is.null(subplot)){
      relations <- apa_config$plot_id_values
    }else{
      relations <- subplot[[agg_column]]
    }
    boundary_length_list <- split(boundary_length_dat, boundary_length_dat[[agg_column]])
    dis_list <- names_to_list(relations)
    if(apa_config$dis_method == "gowdis"){
      dis_fun <- function(...){FD::gowdis(...)}
    }else{
      dis_fun <- apa_config$dis_method
    }
    dis_transform <- apa_config$dis_transform
    ndiv_i <- list()
    for(relation_i in relations){
      ndiv_i[[relation_i]] <- dat_list$tree_dat
      if(inherits(ndiv_i[[relation_i]], "sf")){
        ndiv_i[[relation_i]] <- sf::st_drop_geometry(ndiv_i[[relation_i]])
      }
      ndiv_i[[relation_i]] <- 
        subset(ndiv_i[[relation_i]],
               ndiv_i[[relation_i]][[agg_column]] == relation_i)[, c(agg_column, apa_config$tree_id_column, "apa_size_prop")]
      if(nrow(ndiv_i[[relation_i]]) == 0){
        ndiv_i[[relation_i]][1, ] <- c(list(relation_i), list(NA_character_), list(NA_real_))
      }
      boundary_i <- boundary_length_list[[relation_i]]
      boundary_i_trees <- boundary_i[, paste0(apa_config$tree_id_column, c("", "_bc"))]
      
      trees_present <- 
        stats::na.omit(unique(c(boundary_i_trees[[1]], boundary_i_trees[[2]],
                                ndiv_i[[relation_i]][[2]])))
      trait_val_i <- trait_val[trees_present, , drop = FALSE]
      if(apa_config$scope == "global"){
        trait_val_i <- rbind(trait_val_i, minmax_dummies)
      }
      tree_dis <- as.matrix(dis_transform(dis_fun(trait_val_i)))
      if(all(is.na(boundary_i[[paste0(apa_config$tree_id_column, "_bc")]]))){
        ndiv_i[[relation_i]]$ndiv <- NA_real_
      }else{
        dis_i <- boundary_i
        boundary_i_trees <- subset(boundary_i_trees, !is.na(dis_i[[paste0(apa_config$tree_id_column, "_bc")]]))
        dis_i <- subset(dis_i, !is.na(dis_i[[paste0(apa_config$tree_id_column, "_bc")]]))
        dis_i$.dis <-  tree_dis[as.matrix(boundary_i_trees)]
        dis_i_split <- split(dis_i, dis_i[[apa_config$tree_id_column]])
        ndiv_i_list <- 
          lapply(dis_i_split, 
                 FUN = function(x){
                   data.frame(relation = relation_i,
                              tree_id = x[[apa_config$tree_id_column]][1],
                              ndiv_i = sum(x$boundary_length*x$.dis) / sum(x$boundary_length))
                 })
        ndiv_i_join <- do.call(rbind, ndiv_i_list)
        names(ndiv_i_join) <- c(agg_column, apa_config$tree_id_column, "ndiv")
        match_idx <- match_by(ndiv_i[[relation_i]], ndiv_i_join, by = c(agg_column, apa_config$tree_id_column))
        ndiv_i[[relation_i]]$ndiv <- ndiv_i_join$ndiv[match_idx]
      }
      if(pdiv){
        pdiv_proc <- 
          subset(ndiv_i[[relation_i]],
                 !is.na(ndiv_i[[relation_i]]$apa_size_prop))
        pdiv_tree_id <- pdiv_proc[[apa_config$tree_id_column]]
        tree_dis_pdiv <- 
          tree_dis[pdiv_tree_id, pdiv_tree_id, drop = FALSE]
        pdiv_values <- (pdiv_proc$apa_size_prop %*% tree_dis_pdiv)[1, ]
        ndiv_i[[relation_i]]$pdiv <- ndiv_i[[relation_i]]$apa_size_prop
        if(length(pdiv_values) > 0){
          ndiv_i[[relation_i]]$pdiv[!is.na(ndiv_i[[relation_i]]$apa_size_prop)] <- pdiv_values
        }
      }
      
    }
    ndiv_i <- do.call(rbind, ndiv_i)
    match_idx <- match_by(dat_list$tree_dat, ndiv_i, by = c(agg_column, apa_config$tree_id_column))
    dat_list$tree_dat$ndiv <- ndiv_i[match_idx, ][["ndiv"]]
    dat_list$tree_dat$pdiv <- ndiv_i[match_idx, ][["pdiv"]]
    proc_col <-
      stats::na.omit(c(agg_column, apa_config$agg_class_column, "apa_size",
                       "apa_size_prop", "ndiv", "pdiv"[pdiv]))
    tree_dat_proc <- dat_list$tree_dat[, proc_col]
    if(inherits(tree_dat_proc, "sf")){
      tree_dat_proc <- sf::st_drop_geometry(tree_dat_proc)
    }
    agg_classes <- 
      stats::setNames(c(agg_column, apa_config$agg_class_column),
                      c(names(dat_list)[[1]], apa_config$agg_class_column))
    agg_classes <- stats::na.omit(agg_classes)
    
    dat_list[[1]]$apa_size_prop <- 1
    for(agg_class_i in names(agg_classes)){
      tree_dat_proc_split <- 
        split(tree_dat_proc, paste(tree_dat_proc[[agg_column]], tree_dat_proc[[agg_classes[agg_class_i]]]))
      if(pdiv){
        agg_i_dat <- 
          lapply(tree_dat_proc_split,
                 FUN = function(x){
                   cbind(x[1, unique(c(agg_column, agg_classes[agg_class_i])), drop = FALSE],
                         ndiv = stats::weighted.mean(x$ndiv, x$apa_size, na.rm = TRUE),
                         pdiv = stats::weighted.mean(x$pdiv, x$apa_size, na.rm = TRUE))})
      }else{
        agg_i_dat <- 
          lapply(tree_dat_proc_split,
                 FUN = function(x){
                   cbind(x[1, unique(c(agg_column, agg_classes[agg_class_i])), drop = FALSE],
                         ndiv = stats::weighted.mean(x$ndiv, x$apa_size, na.rm = TRUE))})
      }
      agg_i_dat <- do.call(rbind, agg_i_dat)
      match_columns <- unique(c(agg_column, agg_classes[agg_class_i]))
      match_idx <- match_by(dat_list[[agg_class_i]], agg_i_dat, match_columns)
      dat_list[[agg_class_i]]$ndiv <- agg_i_dat$ndiv[match_idx]
      if(pdiv){
        dat_list[[agg_class_i]]$pdiv <- agg_i_dat$pdiv[match_idx]
      }
    }
    dat_list
  }
