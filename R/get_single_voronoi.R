#' Rasterized version of the 'Weighted Voronoi Diagram'
#'
#' Polygons are rasterized and a Weighted Voronoi Diagram is calculated.
#'
#' @param core Polygon containing a core plot boundary. \code{core} has to
#'   be a length one \code{sfc_POLYGON}-object.
#'
#' @param buffer If available, a polygon containing a buffer plot boundary.
#'  \code{buffer} has to be a length one \code{sfc_POLYGON}-object. Relevant to
#'  decide which assignments of grid cells to points are critical.
#'
#' @param point An \code{sf}-object with a \code{sfc_POINT}geometry column and another
#'   weighting column.
#'
#' @param weight_column Character vector of length one. This column of the
#'   \code{point} dataset will be used as weight. If left out, an unweighted
#'   Voronoi Diagram will be calculated.
#'   
#' @param class_columns Character vector specifying additional columns in
#'  /code{point} that will be added as separate classes to the output
#'  /code{Raster}-object.
#'
#' @param res A numeric vector of length one specifying the resolution of the
#'   output grid
#'   
#' @param warn Issue a warning if the coordinate system is not metric (which is not tested).
#'
#' @return A \code{RasterStack}-object.
#'
#' @keywords internal
get_single_voronoi <- function(core,
                               point,
                               buffer = NULL,
                               weight_column =  NULL,
                               class_columns = NULL,
                               res = 1,
                               warn = TRUE){
  if(is.null(buffer)){
    buffer <- core
  }
  # check input classes
  if(!inherits(core, "sfc_POLYGON")){
    stop("'core' is not an 'sfc_POLYGON'-object.")
  }
  if(!inherits(buffer, "sfc_POLYGON")){
    stop("'buffer' is not an 'sfc_POLYGON'-object.")
  }
  if(!inherits(point, "sf") || !inherits(sf::st_geometry(point), "sfc_POINT")){
    stop("'point' is not an 'sf'-object with a POINT geometry column.")
  }
  
  # check crs
  if(!sf::st_crs(core) == sf::st_crs(buffer)){
    stop("'core' and 'buffer' have to have the same crs.")
  }  
  if(!sf::st_crs(core) == sf::st_crs(point)){
    stop("'core' and 'point' have to have the same crs.")
  }  
  if(warn & !grepl("units=m", sf::st_crs(core)$proj4string)){
    warning("The used crs does not appear to be metric. Algorithms are only tested for metric reference systems.")
  }
  
  # check length of core and buffer
  if(length(core) != 1){
    stop("'core' is not of length one.")
  }
  if(length(buffer) != 1){
    stop("'buffer' is not of length one.")
  }
  
  # check empty features
  if(sf::st_is_empty(core)){
    stop("'core' contains an empty feature")
  }
  if(sf::st_is_empty(buffer)){
    stop("'buffer' contains an empty feature")
  }
  if(any(sf::st_is_empty(point))){
    stop("'point' contains empty features.")
  }
  if(!sf::st_covers(buffer, core, sparse = FALSE)[1, 1]){
    stop("Some part of `core` lies outside of `buffer`. `core` has to be completely covered by `buffer`.")
  }
  # Add weight 1, if no weight specified 
  if(is.null(weight_column)){
    point$.weight <-  1L
    weight_column <- ".weight"
  }
  # Check if weight is numeric
  if(mode(point[[weight_column]]) != "numeric"){
    stop("Weights specified by `", weight_column, "` have to be numeric.")
  }
  # Check for NA in weight and negative weights
  if(any(is.na(point[[weight_column]])) || any(point[[weight_column]] <= 0)){
  }
  
  # check existence of weight and class columns
  if(is.null(point[[weight_column]])){
    stop("'", weight_column, "' is not a column of 'point'.")
  }
  if(!all(class_columns %in% names(point))){
    stop("'", paste(setdiff(class_columns, names(point)), collapse = "', '"),
         "' is not a column of 'point'.")
  }
  raster <- rasterize_polygon(core, res = res)[[1]]
  raster_coords <- get_raster_coordinates(raster)
  
  # The distance of each raster point to the boundary line of the buffer is needed to
  # calculate if assignment of a point to a pixel is critical (affected by an
  # edge effect).
  
  # Convert boundary polygon to linestring
  buffer_coords <- sf::st_coordinates(buffer) 
  # Calculate distance between raster points and boundary line. 
  buffer_dist <- multipoint_linestring_distance(raster_coords, buffer_coords)
  raster_coords <- cbind(raster_coords, buffer_dist)
  
  # Adding weight to points
  point_coords <- cbind(sf::st_coordinates(point), weight =  point[[weight_column]])
  
  # Do the acutal assignment
  point_assignment <- rasterized_weighted_voronoi(raster_coords, point_coords)
  point_assignment <- stats::setNames(point_assignment, c("point_no", "critical"))
  # Select class columns, if specified.
  if(!is.null(class_columns)){
    point_join <- sf::st_set_geometry(point, NULL)
    raster_values <- point_join[point_assignment[[1]], class_columns, drop = FALSE]
    raster_values <- cbind(raster_values, critical = point_assignment[[2]])
  }else{
    raster_values <- as.data.frame(point_assignment)
  }
  # Assign values to raster
  set_raster_values(raster, raster_values)
}
