#' Remove elements from an apa_list
#' 
#' Depending on which function is chosen, either all APA-properties,
#'   all aggregation classes, all subplots or all polygons of APA-patches will
#'   be removed from all datasets in the \code{apa_list}.
#'
#' @template apa_list_arg
#'
#' @return A \code{apa_list} where the respective elements were removed.
#'
#' @example examples/example_apa_drop_x.R
#' 
#' @name apa_drop_x

#' @rdname apa_drop_x
#' @export
#' 
apa_drop_properties <- function(apa_list){
  apa_config <- attr(apa_list, "apa_config")
  drop_properties <- apa_config$apa_properties
  for(data_i in grep("subplot_dat", names(apa_list), value = TRUE, invert = TRUE)){
    select_idx <- !names(apa_list[[data_i]]) %in% drop_properties
    apa_list[[data_i]] <- apa_list[[data_i]][, select_idx , drop = FALSE]
  }
  
  for(subplot_j in names(apa_list$subplot_dat)){
    for(data_i in names(apa_list$subplot_dat[[subplot_j]])){
      select_idx <- !names(apa_list$subplot[[subplot_j]][[data_i]]) %in% drop_properties
      apa_list$subplot_dat[[subplot_j]][[data_i]] <- 
        apa_list$subplot_dat[[subplot_j]][[data_i]][, select_idx , drop = FALSE]
    }
  }
  apa_config_new <- apa_config
  apa_config_new$apa_properties <- NA
  apa_config_new$edge_correction <- NA
  
  apa_config_new$dis_trait_column <- NA
  apa_config_new$dis_method <- NA
  apa_config_new$dis_transform <- NA
  apa_config_new$scope <- NA
  
  apa_config_new <- update_apa_config(apa_list, apa_config_new)
  attr(apa_list, "apa_config") <- apa_config_new
  class(apa_list) <- unique(c("apa_list", class(apa_list)))
  validate_apa_list(apa_list)
}


#' @rdname apa_drop_x
#' @export
apa_drop_agg_class <- function(apa_list){
  apa_config <- attr(apa_list, "apa_config")
    drop_agg_class <- apa_config$agg_class_column
  apa_list$plot_dat$apa_map <-
    lapply(apa_list$plot_dat$apa_map, raster::dropLayer, drop_agg_class)
  apa_list <- apa_list[setdiff(names(apa_list), drop_agg_class)]
  for(subplot_j in names(apa_list$subplot_dat)){
    select_idx <- setdiff(names(apa_list$subplot_dat[[subplot_j]]), drop_agg_class)
    apa_list$subplot_dat[[subplot_j]] <- apa_list$subplot_dat[[subplot_j]][select_idx]
  }
  apa_config_new <- apa_config
  apa_config_new$agg_class_column <- NA
  apa_config_new <- update_apa_config(apa_list, apa_config_new)
  attr(apa_list, "apa_config") <- apa_config_new
  class(apa_list) <- unique(c("apa_list", class(apa_list)))
  validate_apa_list(apa_list)
}

#' @rdname apa_drop_x
#' @export
apa_drop_subplot <- function(apa_list){
  apa_config <- attr(apa_list, "apa_config")
  drop_subplot <- names(apa_config$subplot_id_column)
  apa_list <- apa_list[!names(apa_list) == "subplot_dat"]
  apa_config_new <- apa_config
  apa_config_new$subplot_id_column <- NA
  apa_config_new$radius <- NA
  apa_config_new <- update_apa_config(apa_list, apa_config_new)
  attr(apa_list, "apa_config") <- apa_config_new
  class(apa_list) <- unique(c("apa_list", class(apa_list)))
  validate_apa_list(apa_list)
}

#' @rdname apa_drop_x
#' @export
apa_drop_polygon <- function(apa_list){
  apa_config <- attr(apa_list, "apa_config")
  drop <- function(x){
    x$apa_polygon <- NULL 
    x
  }
  apa_list <- lapply_deep(apa_list, drop)
  apa_config$apa_polygon <- FALSE
  attr(apa_list, "apa_config") <- apa_config
  class(apa_list) <- unique(c("apa_list", class(apa_list)))
  validate_apa_list(apa_list)
}
