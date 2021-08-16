#' Calculate the sizes of APA-patches in an \code{apa_list}
#'
#' For Each dataset in an \code{apa_list}, the sizes of the respective
#' APA-patches will be calculated (absolute and proportional).
#'
#' @template apa_list_arg
#'
#' @inheritParams apa_list
#'
#' @template ref_glatthorn_2021
#'
#' @template ref_gspaltl_2012
#' 
#' @example examples/example_apa_size.R
#'
#' @export
#' 
apa_size <- 
  function(apa_list, 
           edge_correction = "none"){
    apa_config <- attr(apa_list, "apa_config")
    if(!is.na(apa_config$edge_correction) && apa_config$edge_correction != edge_correction){
      stop("There can only be a single edge correction method for one 'apa_list' object.")
    }
    if(edge_correction == "border_tree_exclusion" & !"border_tree" %in% apa_config$apa_properties){
      apa_list <- apa_border_tree(apa_list)
      apa_config <- attr(apa_list, "apa_config")
    }
    if("apa_size" %in% apa_config$apa_properties){
      stop("`apa_size` has already been calculated. Please use `apa_drop_properties()` first to recalculate.")
    }
    apa_config$edge_correction <- edge_correction
    apa_config$apa_properties <- 
      as.vector(stats::na.omit(unique(c(apa_config$apa_properties,
                                        c("apa_size", "apa_size_total", "apa_size_prop", "critical")))))
    apa_map <- apa_list$plot_dat$apa_map
    message("\nAggregating APA-sizes:")
    if(apa_config$core_column == apa_config$buffer_column){
      top_level_plot <- NULL
      top_level_id_column <- NULL
    }else{
      top_level_plot <- apa_list$plot_dat
      sf::st_geometry(top_level_plot) <- apa_config$core_column
      top_level_plot <- top_level_plot[, apa_config$plot_id_column]
      top_level_id_column <- c(plot_dat = apa_config$plot_id_column)
    }
    output <- apa_size_calc(unclass(apa_list),
                            apa_config, apa_map, subplot = top_level_plot,
                            subplot_id_column = top_level_id_column,
                            edge_correction = edge_correction)
    subplot_list <- apa_list$subplot_dat
    for(subplot_i in names(apa_config$subplot_id_column)){
      message("`", subplot_i, "` - aggregating APA-sizes:  ", appendLF = FALSE)
      subplot_dat_i <- apa_list$subplot_dat[[subplot_i]]
      subplot_list[[subplot_i]] <- apa_size_calc(subplot_dat_i, apa_config, apa_map, subplot = subplot_dat_i[[subplot_i]],
                                                 edge_correction = edge_correction, subplot_id_column = apa_config$subplot_id_column[subplot_i])
    }
    message("")
    output$subplot_dat <- subplot_list
    output <- new_apa_list(output, apa_config = apa_config)
    output
  }
apa_size_calc <- 
  function(dat_list, apa_config, apa_map, subplot = NULL, edge_correction = "none", subplot_id_column = NULL){
    if(!is.na(apa_config$edge_correction) && apa_config$edge_correction != edge_correction){
      stop("There can only be a single edge correction method for one 'dat_list' object.")
    }
    if(!edge_correction[1] %in% c("none", "critical", "border_tree_exclusion")){
      stop("`edge_correction` has to be one of `none`, `critical` or `border_tree_exclusion`.")
    }
    if(edge_correction == "critical"){
      critical_layer <- "critical"
    }else{
      critical_layer <- "critical"
    }
    if(edge_correction == "border_tree_exclusion" & !"border_tree" %in% apa_config$apa_properties){
      stop("Border trees have to be identified first with `apa_border_tree()`.")
    }
    if(edge_correction != "border_tree_exclusion" | is.null(subplot)){
      tree_apa_size <- 
        class_area(rst = apa_map,
                   subplot = subplot,
                   layer = apa_config$tree_id_column,
                   critical_layer = critical_layer,
                   plot_id_column = apa_config$plot_id_column,
                   radius = apa_config$radius)
    }else{
      tree_apa_size_values <- 
        class_area(rst = apa_map,
                   subplot = NULL,  # border tree correction uses area size outside the core zone 
                   layer = apa_config$tree_id_column,
                   critical_layer = critical_layer,
                   plot_id_column = apa_config$plot_id_column,
                   radius = apa_config$radius)
      tree_apa_size <- 
        mapply(
          FUN = sf::st_join,
          split(subplot, subplot[[apa_config$plot_id_column]])[apa_config$plot_id_values],
          split(dat_list$tree_dat[, apa_config$tree_id_column], 
                dat_list$tree_dat[[apa_config$plot_id_column]])[apa_config$plot_id_values],
          SIMPLIFY = FALSE)
      tree_apa_size <- do.call(rbind, tree_apa_size)
      tree_apa_size <- sf::st_drop_geometry(tree_apa_size)
      sfc_idx <- sapply(tree_apa_size, inherits, "sfc")
      tree_apa_size <- tree_apa_size[!sfc_idx]
      match_idx <- match_by(tree_apa_size, tree_apa_size_values, by = apa_config$tree_id_column)
      tree_apa_size$area <- tree_apa_size_values$area[match_idx]
    }
    tree_apa_size$apa_size_total <- tree_apa_size$area
    
    names(tree_apa_size)[names(tree_apa_size) == "area"] <- "apa_size"
    if(is.null(subplot)){
      agg_column <- apa_config$plot_id_column
    }else{
      agg_column <- subplot_id_column
    }
    if(edge_correction == "critical"){
      tree_apa_size$apa_size <- tree_apa_size$apa_size * (1 - tree_apa_size$critical)
    }
    if(edge_correction == "border_tree_exclusion"){
      border_tree_id <- dat_list$tree_dat[[apa_config$tree_id_column]][dat_list$tree_dat$border_tree]
      tree_apa_size$border_tree <- tree_apa_size$id_tree_proc %in% border_tree_id
      tree_apa_size$apa_size <- tree_apa_size$apa_size * (1 - tree_apa_size$border_tree)
    }
    apa_size_sum_tmp <- 
      tapply(tree_apa_size$apa_size,
             tree_apa_size[[agg_column]],
             sum)
    apa_size_sum_tmp <- c(apa_size_sum_tmp, recursive = TRUE)
    apa_size_sum_tmp <- apa_size_sum_tmp[tree_apa_size[[agg_column]]]
    tree_apa_size$apa_size_prop <- tree_apa_size$apa_size / apa_size_sum_tmp
    match_columns <- c(agg_column, apa_config$tree_id_column)
    tree_dat_proc <- dat_list$tree_dat
    match_idx <- match_by(tree_dat_proc, tree_apa_size, by = c(agg_column, apa_config$tree_id_column))
    tree_dat_proc$apa_size <- tree_apa_size$apa_size[match_idx]
    tree_dat_proc$apa_size_total <- tree_apa_size$apa_size_total[match_idx]
    tree_dat_proc$apa_size_prop <- tree_apa_size$apa_size_prop[match_idx]
    tree_dat_proc$critical <- tree_apa_size$critical[match_idx]
    tree_dat_proc$critical_area <- tree_dat_proc$critical * tree_dat_proc$apa_size_total
    
    dat_list$tree_dat <- tree_dat_proc
    
    agg_classes <- 
      stats::setNames(c(agg_column, apa_config$agg_class_column),
                      c(names(dat_list)[[1]], apa_config$agg_class_column))
    for(agg_class_i in names(agg_classes)){
      if(is.na(agg_class_i)){
        next
      }
      variables <- c("apa_size", "apa_size_total", "apa_size_prop", "critical_area")
      formula_i <- stats::as.formula(paste0("cbind(", paste(variables, collapse = ", "),
                                            ") ~ ",
                                            agg_column, " + ", agg_classes[agg_class_i]))
      
      agg_i_dat <- stats::aggregate(formula_i, data = tree_dat_proc, 
                                    FUN = function(x)ifelse(all(is.na(x)), NaN, sum(x, na.rm = TRUE)),
                                    na.action = c)
      agg_i_dat$critical <- agg_i_dat$critical_area / agg_i_dat$apa_size_total
      agg_i_dat$critical_area <- NULL
      match_columns <- c(agg_column, agg_classes[agg_class_i])
      match_idx <- match_by(dat_list[[agg_class_i]], agg_i_dat, match_columns)
      variables[4] <- "critical"
      dat_list[[agg_class_i]] <- 
        cbind(dat_list[[agg_class_i]], agg_i_dat[match_idx, variables])
      na_idx <- is.na(dat_list[[agg_class_i]]$apa_size_total)
      #dat_list[[agg_class_i]][na_idx, c("apa_size", "apa_size_total", "apa_size_prop")] <- 0
      dat_list[[agg_class_i]][na_idx, c("apa_size", "apa_size_total")] <- 0
      dat_list[[agg_class_i]][na_idx, "apa_size_prop"] <- NaN
    }
    dat_list
  }
