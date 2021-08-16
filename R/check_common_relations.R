#' @keywords internal
check_common_relations <- function(x, y = NULL, plot_id_column, warn = TRUE){
  if(is.null(y)){
    return(names(x))
  }
  x_name <- as.character(match.call())[2]
  y_name <- as.character(match.call())[3]
  if(inherits(x, 'list')){
    x_label <- names(x)
  } else{
    x_label <- unique(x[[plot_id_column]])
  }
  if(inherits(y, 'list')){
    y_label <- names(y)
  } else{
    y_label <- unique(y[[plot_id_column]])
  }
  common_label <- intersect(x_label, y_label)
  if(length(common_label) == 0){
    stop("There are no common relations in '", x_name, "' and '", y_name,
         "', as specified by 'plot_id_column'.")
  }
  # if(length(x_label) > length(common_label) & warn){
  #   warning("Relations '", paste(setdiff(x_label, common_label), collapse = "', '"),
  #           "' in '", x_name, "' are missing in '", y_name, "' and will be skipped.")
  # }
  # if(length(y_label) > length(common_label) & warn){
  #   warning("Relations '", paste(setdiff(y_label, common_label), collapse = "', '"),
  #           "' in '", y_name, "' are missing in '", x_name, "' and will be skipped.")
  # }
  common_label
}