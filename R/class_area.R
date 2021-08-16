#' @keywords internal
#' 
#' 
class_area <- 
  function(rst, 
           subplot = NULL,
           plot_id_column = NULL,
           layer = NULL,
           radius = 10,
           drop_geometries = TRUE,
           critical_layer = NULL) {
    UseMethod("class_area")
  }

#' @keywords internal
class_area.default <-
  function(rst, 
           subplot = NULL,
           plot_id_column = NULL,
           layer = NULL,
           radius = 10,
           drop_geometries = TRUE,
           critical_layer = NULL){
    
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
    layer <- c(layer, critical_layer)
    if(is.null(layer)){
      layer <- names(rst)
    }else{
      missing_layer <- setdiff(layer, names(rst))
      if(length(missing_layer) > 0){
        stop("'", paste(missing_layer, collapse = "', '", "' is not a layer of 'rst'."))
      }
    }
    rst <- raster::subset(rst, layer, drop = FALSE)
    rst_values <- get_raster_values(rst, return_class = "data.frame")
    rst_values <- lapply(rst_values, as.factor)
    rst_data <- do.call(what = cbind, lapply(rst_values, as.integer)) - 1
    if(is.null(subplot)){
      class_area <- full_ca(x = rst_data)
    }else{
      rst_coord <- get_raster_coordinates(rst)
      rst_data <- cbind(rst_coord, rst_data)
      if(inherits(sf::st_geometry(subplot), "sfc_POLYGON")){
        # passing of spatial objects in R to the cpp boost::geometry functions is done
        # via the wkt format. Not sure if this is the most efficient way to do so but
        # it works. The `digits` argument should be high enough to ensure that
        # coordinates are not truncated. 15 digits should fit for metric CRS.
        subplot_wkt <- sf::st_as_text(sf::st_geometry(subplot), digits = 15) 
        class_area <- subplot_ca(x_wkt = subplot_wkt,  y = rst_data)
        colnames(class_area) <- c(".feature_no", layer, "area")
      }
      if(inherits(sf::st_geometry(subplot), "sfc_POINT")){
        point_coord <- sf::st_coordinates(subplot)
        class_area <- neighborhood_ca(x = point_coord, y = rst_data, radius = radius)
        colnames(class_area) <- c(".feature_no", layer, "area")
      }
    }
    class_area <- as.data.frame(class_area)
    class_area <- subset(class_area, class_area$area > 0)
    class_area$area <- class_area$area * prod(raster::res(rst))
    layer_levels <- lapply(rst_values, levels)
    for(layer_i in layer){
      output_values <- layer_levels[[layer_i]][class_area[, layer_i]]
      if(!raster::is.factor(rst[[layer_i]])){
        output_values <- as.numeric(output_values)
      }
      class_area[[layer_i]] <- output_values
    }
    if(!is.null(critical_layer)){
      agg_formula <- 
        paste0("area ~ .feature_no + ", 
               paste(setdiff(layer, "critical"), collapse = " + "))
      agg_formula <- stats::as.formula(agg_formula)
      class_area_full <- 
        stats::aggregate(agg_formula, data = class_area, FUN =  sum)
      if(all(class_area$critical == 0)){
        class_area <- class_area_full
        class_area$critical <- 0
      }else{
        class_area_crit <- 
          stats::aggregate(agg_formula, data = class_area, sum,
                           subset = class_area$critical == 1)
        names(class_area_crit)[3] <- "critical"
        class_area <- merge(class_area_full, class_area_crit, all = TRUE)
        class_area$critical[is.na(class_area$critical)] <- 0
        class_area$critical <- class_area$critical / class_area$area
      }
    }
    if(is.null(subplot)){
      return(class_area[-1])
    }
    
    output <- as.data.frame(subplot[class_area$.feature_no, ])
    output <- cbind(output, class_area[, -1])
    if(drop_geometries){
      geo_columns <- sapply(output, inherits, "sfc")
      output <- output[!geo_columns]
    }else{
      output <- sf::st_as_sf(output, sf_column_name = attr(subplot, "sf_column"))
    }
    output
  }

#' @keywords internal
class_area.list <- 
  function(rst, 
           subplot = NULL,
           plot_id_column = "raster_name",
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
        class_area(rst_list[[relation_i]],
                   subplot_list[[relation_i]],
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

#' @keywords internal
full_ca <- function(x){
  max_val <- apply(x, 2, max) + 1
  fac_val <- c(1, cumprod(max_val)[-ncol(x)])
  fac_mat <- matrix(rep(fac_val, nrow(x)), ncol = ncol(x), byrow = TRUE)
  count <- tabulate(rowSums(x * fac_mat) + 1)
  req_tab <- cbind(.feature_no = 1, expand.grid(lapply(max_val, seq_len)), area = 0)
  req_tab[1:length(count), "area"] <- count
  req_tab <- req_tab[do.call(what = order, as.data.frame(req_tab)[-ncol(x)-2]), ]
  req_tab
}