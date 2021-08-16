#' Calculate boundary lengths between classes in Raster*-objects
#'
#' Takes a single or a list of Raster*-objects as input and returns border
#' lengths between classes.
#'
#'
#' @keywords internal
boundary_length <- 
  function(rst, 
           layer = NULL,
           subplot = NULL,
           plot_id_column = NULL,
           critical_layer = NULL,
           radius = 10,
           remove_na = FALSE,
           drop_geometries = TRUE,
           border_tree = FALSE) {
    UseMethod("boundary_length")
  }

#'@keywords internal
boundary_length.default <-
  function(rst, 
           subplot = NULL,
           layer = NULL,
           critical_layer = NULL,
           radius = 10,
           remove_na = FALSE,
           drop_geometries = TRUE,
           border_tree = FALSE){
    # check classes
    if(!inherits(rst, c("RasterLayer", "RasterStack", "RasterBrick"))){
      stop("'rst' has to be a 'Raster*'-object or a list of 'Raster*'-objects")
    }
    if(!is.null(subplot) && (!inherits(subplot, "sf") || !inherits(sf::st_geometry(subplot), c("sfc_POINT", "sfc_POLYGON")))){
      stop("'subplot' has to be an 'sf'-object with a 'sfc_POINT'- or 'sfc_POLYGON' geometry.")
    }
    if(!is.null(subplot) && any(sf::st_is_empty(sf::st_geometry(subplot)))){
      stop("No empty geometries allowed in `subplot`.")
    }
    
    rst_res <- raster::res(rst)
    if(!all.equal(rst_res[1], rst_res[2])){
      stop("The resolution of 'res' has to be the same in x and y direction.")
    }
    rst_res <- rst_res[1]
    if(length(layer) > 1){
      stop("'layer' has to point to a single layer.")
    }
    if(is.null(layer)){
      layer <- names(rst)[1]
    }else{
      missing_layer <- setdiff(layer, names(rst))
      if(length(missing_layer) > 0){
        stop("'", paste(missing_layer, collapse = "', '", "' is not a layer of 'rst'."))
      }
    }
    rst_proc <- raster::stack(raster::subset(rst, c(layer, critical_layer)))
    # Extend boundary layer by one pixel and fill them with with `NA`-values. This is
    # necessary to assign edge-pxiels an NA-boundary.
    rst_new_extent <- raster::extend(raster::extent(rst_proc), rep(rst_res, 4))
    rst_proc
    rst_proc <- 
      raster::stack(lapply(raster::unstack(rst_proc), raster::extend, rst_new_extent))
    boundary_cell_info_full <- get_boundary_cell_info(rst_proc[[layer]], radius = radius)
    if(is.null(subplot)){
      boundary_cell_info <- boundary_cell_info_full
    }else{
      boundary_cell_info_full$feature_no <- 0
      boundary_cell_info_subset <- get_boundary_cell_info(rst_proc[[layer]], subset = subplot, radius = radius)
      boundary_cell_info <- rbind(boundary_cell_info_full, boundary_cell_info_subset)
    }
    # two rows are one pair of boundary cells The first one (direction 5) is the
    # focal cell. The second one is the neighbor cell. direction 0 to 3 point
    # out in which direction the boundary is.
    boundary_freq <- 
      count_boundary(boundary_cell_info, rst_proc, layer = layer,
                     critical_layer = critical_layer, remove_na = remove_na,
                     border_tree = border_tree)
    if(!is.null(critical_layer)){
      agg_names <- setdiff(names(boundary_freq), c(critical_layer, "boundary_length"))
      if(all(is.na(boundary_freq[[3]]))){
        boundary_freq_full <- boundary_freq
        boundary_freq_full$critical <- 0
        boundary_freq_full$boundary_length <- 0
      }else{
        boundary_freq_full <- 
          stats::aggregate(stats::as.formula(paste0("boundary_length ~", paste(agg_names, collapse = " + "))),
                           data = boundary_freq,
                           FUN =  sum)
      }
      if(all(is.na(boundary_freq$critical)) || 
         all(is.na(boundary_freq[[3]])) || 
         all(stats::na.omit(boundary_freq$critical) == 0)){
        boundary_freq_crit <- boundary_freq_full
        boundary_freq_crit$critical <- NULL
      }else{
        boundary_freq_crit <- 
          stats::aggregate(stats::as.formula(paste0("boundary_length ~", paste(agg_names, collapse = " + "))),
                           data = boundary_freq, FUN =  sum, subset = boundary_freq[[critical_layer]] == 1)
      }
      names(boundary_freq_crit)[names(boundary_freq_crit) == "boundary_length"] <- "critical"
      boundary_freq <- merge(boundary_freq_full, boundary_freq_crit, all = TRUE)
      boundary_freq$critical[is.na(boundary_freq$critical)] <- 0
      boundary_freq$critical <- boundary_freq$critical / boundary_freq$boundary_length
      boundary_freq$critical[boundary_freq$boundary_length == 0] <- 0 # if there no border on the plot, there is no critical border
    }
    if(is.null(subplot)){
      output <- boundary_freq[, -1]
      return(output)
    }
    boundary_freq <- subset(boundary_freq, boundary_freq$.feature_no != 0)
    output <- as.data.frame(subplot[boundary_freq$.feature_no, ])
    output <- cbind(output, boundary_freq[, -1])
    if(drop_geometries){
      geo_columns <- sapply(output, inherits, "sfc")
      output <- output[!geo_columns]
    }else{
      output <- sf::st_as_sf(output, sf_column_name = attr(subplot, "sf_column"))
    }
    return(output)
  }

#'@keywords internal
boundary_length.list <- 
  function(rst, 
           subplot = NULL,
           plot_id_column = NULL,
           ...){
    if(!is.null(subplot) && (!inherits(subplot, "sf") || !inherits(sf::st_geometry(subplot), c("sfc_POINT", "sfc_POLYGON")))){
      stop("'subplot' has to be an 'sf'-object with a 'sfc_POINT'- or 'sfc_POLYGON' geometry.")
    }
    if(is.null(subplot)){
      rst_list <- check_rst_list(rst)
      common_label <- names(rst_list)
      subplot_list <- names_to_list(common_label, fill = NULL)
    }else{
      common_label <- check_common_relations(rst, subplot, plot_id_column = plot_id_column)
      rst_list <- check_rst_list(rst)[common_label]
      subplot_list <- split(subplot, subplot[[plot_id_column]])[common_label]
    }
    output_list <- names_to_list(common_label)
    for(relation_i in common_label){
      message(relation_i, "  ", appendLF = FALSE)
      output_list[[relation_i]] <- 
        boundary_length(rst = rst_list[[relation_i]],
                        subplot = subplot_list[[relation_i]],
                        ...)
    }
    output <- named_list_to_df(output_list, plot_id_column = plot_id_column)
    if(inherits(output, "sf")){
      # bug in multiple geometry columns: changing sf_column when rbinding
      # and strange bboxes. This is a workaround. TODO: report to edzer pebesma. 
      sf::st_geometry(output) <- attr(subplot, "sf_column") # corrects geometry column
      sfc_idx <- which(sapply(output, inherits, "sfc"))
      for(sfc_i in sfc_idx){
        output[[sfc_i]] <- output[[sfc_i]][1:nrow(output)] # corrects bbox
      }
    }
    output
    
  }
