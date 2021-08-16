#' Efficient transformation of raster classes or patches to polygons
#' @param rst a  \code{Raster*} object
#' @param layer length one character vector. Name of the layer of the
#'   \code{Raster*} object that will be polygonized.
#' @param direction numeric vector (only first element is used). Connectivity of
#'   patches. Either Rook's adjacency (4) or Queen's adjacency (8).
#' @param dissolve_patches length one logical vector. If \code{FALSE},
#'   \code{MULTIPOLYGONS} of class-boundaries are returned. If \code{TRUE},
#'   \code{POLYGONS} of patch-boundaries are returned.
#' @return Ab \code{sf}-object with class labels in the first column and a
#'   \code{geometry} column with class or patch boundaries. if
#'   \code{dissolve_patches} is \code{TRUE}, a column with patch-id's is added.
#' @example examples/example_polygonize_class.R
#' @noRd
polygonize_class <- 
  function(rst,
           layer = NULL,
           direction = c(4, 8),
           plot_id_column = "raster_name",
           geometry_name = "geometry") {
    UseMethod("polygonize_class")
  }

polygonize_class.list <-
  function(rst, 
           layer = NULL,
           plot_id_column = "raster_name",
           ...){
    rst_names <- names(rst)
    if(is.null(rst_names)){
      rst_names <- seq_along(rst_names)
    }
    polygon_list <- names_to_list(rst_names)
    for(rst_i in rst_names){
      message(rst_i, "  ", appendLF = FALSE)
      polygon_list[[rst_i]] <- polygonize_class(rst[[rst_i]], layer, ...)
      polygon_list[[rst_i]][[plot_id_column]] <- rst_i
      polygon_list[[rst_i]] <- polygon_list[[rst_i]][c(3, 1, 2)]
    }
    message("")
    do.call(rbind, polygon_list)
  }

polygonize_class.default <-
  function(rst,
           layer = NULL,
           direction = c(4, 8),
           dissolve_patches = FALSE,
           geometry_name = "geometry"){
    if(!is.null(layer)){
      rst <- raster::subset(rst, layer)
    }else{
      layer <- names(rst)
    }
    if(raster::nlayers(rst) > 1){
      sep <- "__"
      multi_layer <- TRUE
      unique_raster_values <- unique(get_raster_values(rst))
      class_original <- lapply(unique_raster_values, class)
      unique_raster_values <- c(unique_raster_values, recursive = TRUE)
      if(any(grepl(sep, unique_raster_values)) | any(grepl(sep, layer))){
        stop("`__` (double underscore) is not allowed to occur in layer-names or -values`")
      }
      rst <- raster::setValues(x = raster::raster(rst), 
                               values = factor(apply(raster::values(rst), MARGIN = 1, 
                                                     FUN = paste, collapse = sep)))
      layer_original <- layer
      layer <- paste(layer, collapse = sep)
      names(rst) <- layer
      #stop("Please provide either a `RasterLayer`-object or specify which (individual) layer to use in `layer`.")
    }else{
      multi_layer <- FALSE
    }
    # Extend boundary layer by one pixel and fill them with with `NA`-values. This is
    # necessary to assign edge-pxiels an NA-boundary.
    rst_res <- raster::res(rst)
    if(!all.equal(rst_res[1], rst_res[2])){
      stop("The resolution of 'res' has to be the same in x and y direction.")
    }
    rst_res <- rst_res[1]
    rst_new_extent <- raster::extend(raster::extent(rst), rep(rst_res, 4))
    rst_proc <- raster::extend(rst, rst_new_extent)
    
    rst_boundary <- raster::as.matrix(rst_proc)
    
    boundary_mat <- get_boundary_pixels(x = rst_boundary)
    boundary_cell_info <- get_boundary_cell_info(rst = rst_proc)
    
    boundary_cell_no <- 
      boundary_cell_info$col_n + (boundary_cell_info$row_n - 1) * ncol(rst_proc)
    
    # rownumber of `layer_values` corresponds to cell_no (because of `remove_na = FALSE`).
    layer_values <- get_raster_values(rst_proc, remove_na = FALSE)
    boundary_dat <- 
      cbind(class = layer_values[boundary_cell_no, , drop = FALSE],
            direction = boundary_cell_info$direction,
            row = boundary_cell_info$row_n,
            col = boundary_cell_info$col_n)
    boundary_center_sub <- subset(boundary_dat, direction == 5)
    boundary_neighbor_sub <- subset(boundary_dat, direction != 5)
    boundary_wide_dat <- 
      cbind(boundary_center_sub[c(layer, "row", "col")],
            boundary_neighbor_sub["direction"])
    boundary_wide_dat <- subset(boundary_wide_dat, !is.na(boundary_wide_dat[[layer]]))
    boundary_wide_dat$direction_t <- c(1, 2, 3, 0)[boundary_wide_dat$direction + 1]
    boundary_wide_dat$row1 <- boundary_wide_dat$row * 2 + c(-1, -1, +1, +1)[boundary_wide_dat$direction + 1]
    boundary_wide_dat$col1 <- boundary_wide_dat$col * 2 + c(-1, +1, +1, -1)[boundary_wide_dat$direction + 1]
    boundary_wide_dat$row2 <- boundary_wide_dat$row * 2 + c(-1, +1, +1, -1)[boundary_wide_dat$direction + 1]
    boundary_wide_dat$col2 <- boundary_wide_dat$col * 2 + c(+1, +1, -1, -1)[boundary_wide_dat$direction + 1]
    boundary_wide_dat <- 
      boundary_wide_dat[c(layer, "direction_t", "row1", "col1", "row2", "col2")]
    sort_idx <- do.call(order, boundary_wide_dat[, c(layer), drop = FALSE])
    boundary_wide_dat <- boundary_wide_dat[sort_idx, ]
    boundary_list <- 
      split(boundary_wide_dat[, -1],boundary_wide_dat[[layer]])
    boundary_list <- 
      lapply(boundary_list, as.matrix)
    # The actual polygonization
    boundary_poly_coord <- 
      lapply(boundary_list, polygonize_class_cpp, direction = direction[1])
    
    #Transformation of the coordinates
    rst_ymax <- raster::ymax(rst_proc) - rst_res / 2
    rst_xmin <- raster::xmin(rst_proc) + rst_res / 2
    boundary_poly_coord <-  
      lapply(boundary_poly_coord, lapply, lapply, 
             function(x)cbind(rst_xmin + rst_res * ((x[, 2]) / 2 - 1),
                              rst_ymax - rst_res * ((x[, 1]) / 2 - 1)))
    
    # Transformation of coordinates to multipolygons (sfc) and then to an sf-object
    boundary_multipoly <- 
      lapply(boundary_poly_coord, sf::st_multipolygon)
    boundary_sfc <- sf::st_sfc(boundary_multipoly, crs = raster::crs(rst, asText = TRUE))
    boundary_sfc <-lwgeom::lwgeom_make_valid(boundary_sfc)
    
    boundary_sf <- unique(boundary_wide_dat[, layer, drop = FALSE])
    boundary_sf[[geometry_name]] <- boundary_sfc
    sf::st_geometry(boundary_sf) <- geometry_name
    sf::st_agr(boundary_sf) <- "constant"
    boundary_sf <- sf::st_as_sf(boundary_sf)
    if(multi_layer){
      boundary_val_original <- do.call(rbind, strsplit(boundary_sf[[layer]], sep))
      boundary_val_original <- as.data.frame(boundary_val_original)
      names(boundary_val_original) <- layer_original
      boundary_val_original <- 
        as.data.frame(mapply(boundary_val_original, class_original, FUN = methods::as))
      boundary_sf <- 
        sf::st_as_sf(cbind(boundary_val_original, boundary_sf[, "geometry"]))
      
    }
    boundary_sf
  }
