#' @keywords internal
border_class <- 
  function(rst, layer) {
    UseMethod("border_class")
  }

border_class.default <- function(rst, layer = NULL){
  if(is.null(layer)){
    layer <- names(rst[[1]])
  }
  rst <- raster::subset(rst, layer)
  rst_res <- raster::res(rst)
  if(!all.equal(rst_res[1], rst_res[2])){
    stop("Resolution has to be identical in `x` and `y` direction.")
  }
  rst_res <- rst_res[1]
  rst_new_extent <- raster::extend(raster::extent(rst), rep(rst_res, 4))
  rst <- raster::extend(rst, rst_new_extent)
  rst_mat <- raster::as.matrix(rst)
  na_idx <- is.na(rst_mat)
  na_col <- col(rst_mat)[na_idx]
  na_row <- row(rst_mat)[na_idx]
  na_rowcol <- cbind(na_row, na_col)
  # queens case:
  offset_list <-
    list(c(-1, -1), c(0, -1), c(1, -1), c(-1, 0), c(1, 0), c(-1, 1), c(0, 1), c(1, 1))
  # rooks case
  offset_list <- offset_list[sapply(offset_list, `%in%`, x = 0)]
  
  border_cell_list <- lapply(offset_list, function(x) na_rowcol + rep(x, each = nrow(na_rowcol)))
  border_cells <- do.call(rbind, border_cell_list)
  inner_idx <- !border_cells[, 1] %in% c(0, nrow(rst_mat) + 1) &
    !border_cells[, 2] %in% c(0, ncol(rst_mat) + 1)
  border_cells <- border_cells[inner_idx, ]
  border_cell_mat <- matrix(rep(FALSE, length(rst_mat)), ncol = ncol(rst_mat), nrow = nrow(rst_mat))
  border_cell_mat[border_cells] <- TRUE
  border_cell_mat[is.na(rst_mat)] <- FALSE
  rst_values <- get_raster_values(rst, remove_na = FALSE)
  border_class <- unique(rst_values[t(border_cell_mat), ])
  rst_values <- unique(subset(rst_values, !is.na(rst_values[[1]])))
  rst_values$border_class <- rst_values[[1]] %in% border_class
  rownames(rst_values) <- rst_values[[1]]
  rst_values
}

border_class.list <- 
  function(rst, ...){
    rst_names <- names(rst)
    if(is.null(rst_names)){
      rst_names <- seq_along(rst_names)
    }
    output_list <- names_to_list(rst_names)
    message("Looking up border classes:")
    for(rst_i in rst_names){
      message(rst_i, "  ", appendLF = FALSE)
      output_list[[rst_i]] <- border_class(rst[[rst_i]], ...)
    }
    message("")
    output <- do.call(rbind, output_list)
    rownames(output) <- output[[1]]
    output
  }
