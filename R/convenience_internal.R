#' @keywords internal 
names_to_list <- function(x, fill = NA, output_class = "list"){
  output <- structure(rep(list(fill), length(x)), names = x)
  if(output_class == "data.frame"){
    output <- as.data.frame(output)
  }
  output
}

#' @keywords internal 
named_list_to_df <- 
  function(x, plot_id_column){
    x_freq <- sapply(x, nrow)
    x_label <- rep(names(x_freq), times = x_freq)
    output <- do.call(rbind, x)
    if(!is.null(output[[plot_id_column]])){
      if(all.equal(as.character(output[[plot_id_column]]), x_label)){
        return(output)
      }
      stop("Names of 'x' and 'x[[plot_id_column]]' are not identical.")
    }
    x_label_df <- stats::setNames(data.frame(x_label, stringsAsFactors = FALSE), plot_id_column)
    output <- cbind(x_label_df, output)
    output
  }


#' @keywords internal 
match_by <- function(to, from, by){
  if(inherits(to, "sf")){
    to <- sf::st_drop_geometry(to)
  }
  if(inherits(from, "sf")){
    from <- sf::st_drop_geometry(from)
  }
  match_cat_to <- do.call(paste, to[, by, drop = FALSE])
  match_cat_from <- do.call(paste, from[, by, drop = FALSE])
  match_idx <- match(match_cat_to, match_cat_from)
}

#' @keywords internal 
lapply_deep <- function(.x, .f, .p = inherits, .what = "data.frame", .first = identity, .other = identity, ...){
  if(ifelse(all.equal(.p, inherits), .p(.x, .what), .p(.x))){
    return(.f(.x, ...))
  }
  if(!is.null(.first)){
    .x <- .first(.x)
    .first <- NULL
  }else{
    .x <- .other(.x)
  }
  lapply(.x, lapply_deep, .f = .f, .p = .p, .what = .what, .first = NULL, .other = .other, ...)
}

#' @keywords internal 
mapply_deep <- function(.x, .f, .p = inherits, .what = "data.frame", ...){
  if(ifelse(all.equal(.p, inherits), .p(.x[[1]], .what), .p(.x[[1]]))){
    return(.f(.x, ...))
  }else{
    output <- list()
    for(i in seq_along(.x[[1]])){
      .x_i <- lapply(.x, '[[', i)
      output[[names(.x[[1]])[i]]] <- 
        mapply_deep(.x_i, .f = .f, .p = .p, .what = .what, ...)
    }
  }
  output
}

#' @keywords internal 
combine_df_cols <- function(x, value_columns){
  if(nrow(x[[1]]) == 0){
    return(data.frame())
  }
  output_list <- lapply(x, function(x){x[setdiff(names(x), value_columns)]})
  check_id <- all(sapply(output_list[-1], all.equal, target = output_list[[1]]))
  if(!check_id){
    stop("Cannot combine becuase id columns are not identical")
  }
  output <- output_list[[1]]
  for(column_j in value_columns){
    output_column_j <- lapply(x, `[[`, column_j)
    output_column_j <- do.call(cbind, output_column_j)
    output_column_j <- unname(split(output_column_j, row(output_column_j)))
    output[[column_j]] <- output_column_j
  }
  output
}

#' @keywords internal 
drop_columns_except <- function(x, keep_columns){
  if(inherits(x, "sf")){
    x <- sf::st_drop_geometry(x)
  }
  x <- x[, intersect(names(x), keep_columns), drop = FALSE]
}

#' @keywords internal 
add_col <- function(x, col){
  for(col_i in col){
    if(length(x[[2]][[col_i]]) > 0){
      x[[1]][[col_i]] <- x[[2]][[col_i]]
    }
  }
  x[[1]]
}


#' @keywords internal 
shorten_string <- 
  function(x, head_n, shorten_to){
    x_string <- utils::head(x, head_n)
    x_append <- ifelse(nchar(x_string) > shorten_to, "...", rep("", length(x_string)))
    x_string <- paste0(substr(x_string, 1, shorten_to), x_append)
    x_string
  }
