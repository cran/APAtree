# Returns coordinates of a raster* object as matrix. NA values can optionally be
# removed.
#' @keywords internal
get_raster_coordinates <- function(raster, remove_na = TRUE, as_sfc = FALSE){
  
  raster_coordinates <- 
    raster::xyFromCell(raster, cell = 1:raster::ncell(raster))
  
  if(remove_na){
    raster_coordinates <- 
      raster_coordinates[!is.na(raster::getValues(raster[[1]])), , drop = FALSE]
  }
  
  if(!as_sfc){
    return(raster_coordinates)
  }
  
  crs <- sf::st_crs(raster::crs(raster))
  
  sf::st_geometry(sf::st_as_sf(as.data.frame(raster_coordinates),
                               coords = 1:2, 
                               crs = crs))
}

#' @keywords internal
get_raster_values <- function(rst,
                              layer = NULL,
                              remove_na = TRUE,
                              return_class = c("data.frame", "matrix"),
                              raster_factor_as = "character",
                              raster_logical_as = "integer"){
  if(is.null(layer)){
    layer <- names(rst)
  }
  rst <- raster::subset(rst, layer)
  rst_values <- as.matrix(raster::getValues(rst))
  colnames(rst_values) <- layer
  if(remove_na){
    na_values <- is.na(raster::getValues(rst[[1]]))
    rst_values <- rst_values[!na_values, , drop = FALSE]
  }
  if(raster_logical_as == "integer"){
    mode(rst_values) <- "numeric"
  }
  if(return_class[1] == "matrix"){
    return(rst_values)
  }
  if(raster_factor_as == "integer"){
    return(as.data.frame(rst_values))
  }
  raster_levels <- raster::levels(rst)
  # unnesting of strange raster::levels output if variables of different types are mixed
  nested_idx <- sapply(raster_levels, length) == 1
  raster_levels[nested_idx] <- 
    lapply(raster_levels[nested_idx], `[[`, 1)
  if(is.null(raster_levels)){
    raster_levels <- rep(list(list()), length(layer))
  }
  names(raster_levels)  <-  layer
  output_list <- structure(vector('list', length(layer)), names = layer)
  for(layer_i in layer){
    layer_values <- rst_values[, layer_i]
    if(length(raster_levels[[layer_i]]) > 0){
      layer_values <- raster_levels[[layer_i]][layer_values, 2]
    }
    output_list[[layer_i]] <- layer_values
  }
  if(raster_factor_as == "character"){
    fac_idx <- sapply(output_list, is.factor)
    output_list[fac_idx] <- lapply(output_list[fac_idx], as.character)
  }
  as.data.frame(output_list, stringsAsFactors = FALSE)
}



#' @keywords internal
rst_reclassify_chr <- function(rst, reclass_df){
  if("critical" %in% names(rst)){
    rst_critical <- rst[["critical"]]
    rst <- raster::subset(rst, setdiff(names(rst), "critical"))
  }else{
    rst_critical <- NULL
  }
  rst_values <- get_raster_values(rst[[names(reclass_df)[1]]], remove_na = FALSE)
  match_idx <- match(rst_values[[1]], reclass_df[[1]])
  rst_values_classified <- reclass_df[match_idx, -1, drop = FALSE]
  rst_values_classified <- lapply(rst_values_classified, factor)
  rst_reclassified <- 
    lapply(rst_values_classified, raster::setValues, x = raster::raster(rst))
  rst_list <- stats::setNames(raster::unstack(raster::stack(rst)), names(rst))
  rst_reclassified <- raster::stack(c(rst_list, rst_reclassified, critical = rst_critical))
}


#' Replacement of non-NA values of a raster
#'
#' Takes \code{raster*} objects as input and replaces all non-NA values with
#' values specified in a \code{data.frame}. Multiple columns will return
#' \code{RasterStack}'s.
#'
#' @param raster A \code{raster*} object. All non-NA values of \code{raster}
#'   will be replaced. If \code{raster} is a \code{RasterStack}, each layer has
#'   to have identical \code{na}-cells.
#'
#' @param data A \code{data.frame} with as many rows as there are non-NA values
#'   in \code{raster}.
#'
#' @return A \code{RasterLayer} if \code{data} has only one column. A
#'   \code{RasterBrick} if multiple columns were specified in \code{data}.
#'   Character and factor data will lead to the creation of factor layers.
#' @keywords internal
set_raster_values <- function(raster, data){
  # remark: raster* objects are always double matrices. There are no integer matrices. 
  na_values <- is.na(raster::getValues(raster[[1]]))
  data_fac <- sapply(data, is.factor)
  data[, data_fac] <- lapply(data[, data_fac, drop = FALSE], as.character)
  
  data_with_na <- 
    lapply(data, function(x){replace(rep(NA, length(na_values)), !na_values, x)})
  data_with_na <- do.call(what = data.frame, data_with_na)
  
  data_with_na_char <- sapply(data_with_na, is.character)
  data_with_na[, data_with_na_char] <-
    lapply(data_with_na[, data_with_na_char, drop = FALSE], as.factor)
  
  dummy_raster <- raster::raster(raster)
  raster_layer_list <- 
    lapply(data_with_na, function(x){raster::setValues(dummy_raster, x)})
  
  if(length(raster_layer_list) == 1){
    return(raster_layer_list[[1]])
  }
  raster::stack(raster_layer_list)
}


