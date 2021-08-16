#' @keywords internal
get_boundary_cell_info <- function(rst, subset = NULL, radius = NULL){
  rst_mat <- raster::as.matrix(rst)
  if(typeof(rst) == "double"){
    stop("The 'boundary_length'-algorithm only works for categorical data. The selected 'layer' is of type double.")
  }
  
  # Extract boundary pixels (including edge- (NA-)boundaries)
  full_boundary_cell_info <- get_boundary_pixels(x = rst_mat)
  colnames(full_boundary_cell_info) <- c("direction", "row_n", "col_n")
  # two rows are one pair of boundary cells The first one (direction 5) is the
  # focal cell. The second one is the neighbor cell. direction 0 to 3 point
  # out in which direction the boundary is.
  if(is.null(subset)){
    boundary_cell_info <- as.data.frame(cbind(feature_no = 1, full_boundary_cell_info))
    return(boundary_cell_info)
  }
  full_boundary_cell_coordinates <-
    cbind(
      raster::xFromCol(rst, full_boundary_cell_info[, "col_n"]),
      raster::yFromRow(rst, full_boundary_cell_info[, "row_n"]))
  full_boundary_cell_info <- 
    cbind(full_boundary_cell_coordinates, full_boundary_cell_info)
  if(inherits(sf::st_geometry(subset), "sfc_POINT")){
    subset_coord <- sf::st_coordinates(subset)
    boundary_cell_info <-
      as.data.frame(get_neighborhood_boundary_pixels(subset_coord, full_boundary_cell_info, radius = radius))
    names(boundary_cell_info) <- c("feature_no", "direction", "row_n", "col_n")
  }
  if(inherits(sf::st_geometry(subset), "sfc_POLYGON")){
    subset_wkt <- sf::st_as_text(sf::st_geometry(subset), digits = 15)
    boundary_cell_info <-
      as.data.frame(get_subplot_boundary_pixels(subset_wkt, full_boundary_cell_info))
    names(boundary_cell_info) <- c("feature_no", "direction", "row_n", "col_n")
  }
  boundary_cell_info
}
  
  