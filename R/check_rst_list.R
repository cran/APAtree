#' @keywords internal
check_rst_list <- function(rst){
  rst_idx <- sapply(rst, inherits, c("RasterLayer", "RasterStack", "RasterBrick"))
  if(!any(rst_idx)){
    stop("There are no 'Raster*'-objects in 'rst'.")
  }
  if(!all(rst_idx)){
    warning("Some objects in 'rst' are not of one of the 'Raster*'-classes.
              They will not be processed.")
  }
  rst_label <- names(rst)
  if(is.null(rst_label)){
    rst_label <- paste0("rst_", 1:length(rst_label))
  }
  rst[rst_idx]
}