#' Assess which trees in APA-maps are border trees.
#'
#' Border trees are trees whose APA-patches touch the plot-border.
#'
#' @template apa_list_arg
#'
#' @return An \code{apa_list} where \code{apa_list$tree_dat} will have an
#'   additional column that specifies which trees are border trees.
#'
#' @export
#' @example examples/example_border_tree.R
apa_border_tree <- function(apa_list){
  apa_config <- attr(apa_list, "apa_config")
  if(apa_config$core_column != apa_config$buffer_column){
    stop("To assign border trees, 'core' and 'buffer' polygons have to be identical.")
  }
  message("\nIdentifying border trees.")
  border_tree <- border_class(apa_list$plot_dat$apa_map)
  apa_list$tree_dat$border_tree <- 
    border_tree[apa_list$tree_dat[[apa_config$tree_id_column]], 2]
  for(subplot_i in names(apa_config$subplot_id_column)){
    apa_list$subplot_dat[[subplot_i]]$tree_dat$border_tree <- 
      border_tree[apa_list$subplot_dat[[subplot_i]]$tree_dat[[apa_config$tree_id_column]], 2]
  }
  apa_config$apa_properties <- 
    as.vector(stats::na.omit(unique(c(apa_config$apa_properties, "border_tree"))))
  apa_list <- new_apa_list(apa_list, apa_config = apa_config)
  apa_list
}