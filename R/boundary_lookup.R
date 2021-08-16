#' @keywords internal
count_boundary <- function(boundary_cell_info, rst, layer, critical_layer, remove_na, border_tree){
  rst_res <- raster::res(rst)
  if(!all.equal(rst_res[1], rst_res[2])){
    stop("The resolution of 'res' has to be the same in x and y direction.")
  }
  rst_res <- rst_res[1]
  
  boundary_cell_no <- 
    boundary_cell_info$col_n + (boundary_cell_info$row_n - 1) * ncol(rst)
  
  # rownumber of `layer_values` corresponds to cell_no (because of `remove_na = FALSE`).
  layer_values <- get_raster_values(rst, return_class = "data.frame", remove_na = FALSE)
  boundary_dat <- 
    cbind(.feature_no = boundary_cell_info$feature_no,
          direction = boundary_cell_info$direction,
          layer_values[boundary_cell_no, , drop = FALSE])
  # data.frame with one row per boundary element:
  boundary_center_sub <- subset(boundary_dat, direction == 5)[, c(".feature_no", layer, critical_layer), drop = FALSE]
  boundary_neighbor_sub <- subset(boundary_dat, direction != 5)[, c(layer, critical_layer), drop = FALSE]
  # When a `critical_layer` is selected, a single categorie for boundaries with
  # no critical cells is created.
  if(!is.null(critical_layer)){
    crit_boundary <- 
      as.numeric(boundary_center_sub[[critical_layer]] | boundary_neighbor_sub[[critical_layer]])
    boundary_center_sub <- boundary_center_sub[, c(".feature_no", layer), drop = FALSE]
    boundary_neighbor_sub <- boundary_neighbor_sub[, layer, drop = FALSE]
  }else{
    crit_boundary <- NULL
  }
  names(boundary_center_sub)[2] <- layer
  names(boundary_neighbor_sub)[1] <- paste0(layer, "_bc")
  boundary_wide_dat <- 
    cbind(boundary_center_sub, 
          boundary_neighbor_sub)
  boundary_wide_dat$critical <- crit_boundary
  boundary_wide_dat <- subset(boundary_wide_dat, !is.na(boundary_wide_dat[, 2]))
  if(border_tree){
    border_tree_label <- unique(boundary_wide_dat[is.na(boundary_wide_dat[, 3]), 2])
  }
  # removal of NA-boundaries (if selected via argument `na_remove = TRUE`):
  if(remove_na){
    boundary_wide_dat <- 
      subset(boundary_wide_dat, !is.na(boundary_wide_dat[, 3]))
  }
  boundary_freq <- 
    as.data.frame(table(boundary_wide_dat, useNA = "ifany"),
                  responseName = "boundary_length",
                  stringsAsFactors = FALSE)
  if(border_tree){
    boundary_freq$border_tree = boundary_freq[, 2] %in% border_tree_label
  }
  boundary_freq$.feature_no <- as.numeric(boundary_freq$.feature_no)
  boundary_freq <- subset(boundary_freq, boundary_length > 0)
  if(!is.null(critical_layer)){
    boundary_freq$critical <- as.numeric(boundary_freq$critical)
  }
  if(nrow(boundary_freq) > 0){
    boundary_freq <- boundary_freq[do.call(order, boundary_freq[, -4]), ]
  }
  boundary_freq$boundary_length <- boundary_freq$boundary_length * rst_res
  boundary_freq
}
