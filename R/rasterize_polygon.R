#' Rasterize one or multiple polygons. Takes sf or sp polygons as inputs and
#' returns a list with one \code{RasterLayer} per polygon.
#' @keywords internal
rasterize_polygon <- 
  function(polygon,
           res = 1, 
           labels = "layer"){
  geo_empty_idx <- sf::st_is_empty(polygon)
  
  polygon_crs <- sf::st_crs(polygon)$proj4string
  
  non_empty_geo = sf::st_geometry(polygon)[!geo_empty_idx]
  polygon_sp <- lapply(non_empty_geo, methods::as, "Spatial")
  empty_raster <- lapply(polygon_sp, raster::raster, res = res)
  empty_raster <- lapply(empty_raster, function(x){raster::crs(x) <- polygon_crs; x})
  raster_list <- 
    mapply(FUN = raster::rasterize,
           polygon_sp,
           empty_raster,
           getCover = FALSE)
  raster_list <- ifelse(geo_empty_idx, NA, raster_list)
  if(inherits(polygon, "sf")){
    if(length(labels) == 1){
      poly_labels <-  polygon[[labels]]
    }else if(is.null(labels)){
      poly_labels <- paste0("p", seq_len(nrow(polygon)))
    }else{
      poly_labels = labels
    }
  } else {
    poly_labels = labels
  }
  stats::setNames(raster_list, poly_labels)
}
