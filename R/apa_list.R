#' Calculating APA-maps
#'
#' Calculate maps of the area potentially available to trees (APA-maps) out of
#' tree inventory data. All vector data that is used or provided by the
#' \code{APAtree}-package is stored in \code{data.frame}s of the additional
#' \code{sf} class of the [`sf`][sf::sf]-package (referred to as
#' \code{sf-data.frame}). APA-maps are raster data that are stored in objects of
#' the [`raster`][raster::raster-package]-package.
#'
#' APA-maps are derived by calculating a rasterized version of weighted
#' Voronoi-Diagrams of tree coordinates. To each tree a weight may be assigned
#' that reflects the competitive ability of the tree.
#'
#' @param plot_dat An \code{sf-data.frame} with plot-level data. \code{plot.dat}
#'   must contain a column that specifies the unique id of the plot. the column
#'   name of the plot id is specified with the \code{plot_id_column}-argument.
#'   The plot id is used to relate plot- and tree-data. A column with
#'   \code{POLYGON}-data that specify the outline of the plots has to be
#'   provided.
#'
#' @param tree_dat An \code{sf}-\code{data.frame} with tree data.
#'   \code{tree.dat} must contain a unique id of each tree. The column name of
#'   the tree id is specified with \code{tree_id_column}. Plot and tree-data are
#'   related by the id of the plot, which has to be a column in both datasets.
#'   The geometry type of \code{tree_dat} has to be \code{POINT}-data that
#'   specifies tree coordinates.
#'
#' @param plot_id_column Column name of \code{plot_dat} specifying unique plot
#'   id's.
#'
#' @param tree_id_column Column name of \code{tree_dat} specifying unique tree
#'   id's.
#'
#' @param weight_column Column name of \code{tree_dat} that specifies a variable
#'   to be used as weighting factor for the calculation of APA-maps.
#'
#' @param agg_class_column An optional vector of column names of \code{tree_dat}
#'   that specify variables that are used as additional grouping variable (for
#'   example tree species) to aggregate apa-properties (see [apa_add_agg_class]
#'   for more details).
#'
#' @param core_column Column name of \code{plot_dat} specifying a
#'   [`sfc`][sf::sfc]-column with \code{POLYGON}-data of plot boundaries.
#'   Default is the active geometry of \code{plot_dat}.
#'
#' @param buffer_column Column name of \code{plot_dat} specifying a
#'   [`sfc`][sf::sfc]-column with \code{POLYGON}s of boundaries of a buffer-zone
#'   around the core plot. To specify a buffer is only necessary if trees
#'   outside of the core area of a plot were sampled (plus-sampling) and the
#'   \code{edge-correction} method is \code{"critical"}. Defaults to
#'   \code{core_column}.
#'
#' @param res single number specifying the resolution of the APA-maps (identical
#'   in x- and y-direction).
#'
#' @param subplot_dat A named list of \code{sf-data.frames} with additional data
#'   of subplots. If the active geometry column of \code{plot_dat} is
#'   \code{POLYGON}-data, calculation of apa-properties will be done for the
#'   area within these \code{POLYGON}s. For \code{POINT}-data, circular subplots
#'   are created around the point coordinates.
#'
#' @param subplot_id_column A named vector specifying the id columns of
#'   \code{subplot_dat}.
#'
#' @param radius The radius of circular subplots that should be used if a
#'   dataset in in \code{subplot_dat} contains \code{POINT}-data.
#'
#' @param apa_properties A vector specifying the apa-properties that will be
#'   calculated. May be any combination of \code{"border_tree"},
#'   \code{"apa_size"}, \code{"ndiv"} and \code{"pdiv"}.
#'
#' @param edge_correction which of the implemented edge correction method should
#'   be applied when calculating apa-properties (\code{"none"},
#'   \code{"critical"} or \code{"border_tree_exclusion"}.
#'
#' @param dis_trait_column A list containing combinations of traits that will
#'   be used to estimate dissimilarity between trees when calculating
#'   apa-properties. Refers to column names of \code{tree_dat}.
#'
#' @param dis_method Which method should be used to estimate dissimilarity
#'   between trees. If \code{dis_method} is \code{gowdis}, [FD::gowdis()] will
#'   be used to calculate Gower dissimilarity. Alternatively, \code{dis_method}
#'   may be any \code{function} that calculates a dissimilarity matrix out of a
#'   \code{data.frame} with trait values.
#'
#' @param dis_transform A \code{function} to transform dissimilarities. Defaults
#'   to \code{sqrt}.
#'
#' @param scope Should scaling of the dissimilarity be done at \code{"global"}
#'   or at \code{"local"} level? If \code{dis_method} scales dissimilarity
#'   between trees according to the range of occurring values in the dataset (as
#'   done by \code{"gowdis"}), a \code{"global"} \code{scope} will use the range
#'   of all values in \code{tree_dat}. Any other \code{scope} will use the range
#'   of values at plot-level to scale dissimilarity.
#'
#' @param apa_polygon logical, specifies if \code{POLYGON}s of APA-patches be
#'   should be added to the datasets. Defaults to \code{TRUE}.
#'
#' @return \code{apa_list} returns an object of class \code{"apa_list"}. An
#'   \code{apa_list} is a \code{list} of \code{data.frame}s with at least two
#'   elements:
#'
#'   * \code{apa_list$plot_dat}: A \code{data.frame} with the original
#'   plot-level data is stored with an additional column \code{apa_map} that is
#'   a list containing all APA-maps. APA-maps of individual plots are
#'   represented by \code{RasterStack} objects. If additional APA-properties
#'   are calculated at stand-level, they are appended to this \code{data.frame}
#'   as well.
#'
#'   * \code{apa_list$tree_dat}: A \code{data.frame} with the original
#'   tree-level data Additional tree-level APA-characteristics are added to
#'   \code{tree_dat} as separate columns.
#'
#'   * aggregation classes (optional): if one or multiple aggregation classes
#'   (e.g., tree species) were added to the \code{apa_list}, additional
#'   \code{data.frame}s at class-level will be added..
#'
#'   * \code{apa_list$subplot_dat} (optional): If subplots were specified, all
#'   data at subplot-level are stored here.
#' 
#' @template ref_glatthorn_2021
#' @template ref_roemisch_1995
#' @template ref_gspaltl_2012
#' 
#' @example examples/example_apa_list.R
#'
#' @export
#' @md
#' 
apa_list <-
  function(plot_dat,
           tree_dat,
           plot_id_column,
           tree_id_column,
           weight_column,
           agg_class_column = NULL,
           core_column = attr(plot_dat, "sf_column"),
           buffer_column = core_column,
           res = 1,
           subplot_dat = NULL,
           subplot_id_column = NULL,
           radius = 10,
           apa_properties = NA,
           edge_correction = "none",
           dis_trait_column = NULL,
           dis_method = "gowdis",
           dis_transform = sqrt,
           scope = "global",
           apa_polygon = TRUE){
    if(any(duplicated(c(tree_id_column, plot_id_column, subplot_id_column)))){
      stop("`plot_id_column` and id columns of `tree_dat` or `subplot_dat` must not be identical.")
    }
    if(is.null(tree_dat[[tree_id_column]])){
      stop("`tree_dat[[tree_id_column]]` does not exist.")
    }
    if(mode(tree_dat[[tree_id_column]]) != "character"){
      stop("`tree_dat[[tree_id_column]]` has to be a `character` vector.")
    }
    if(is.null(tree_dat[[plot_id_column]])){
      stop("`tree_dat[[plot_id_column]]` does not exist.")
    }
    if(mode(tree_dat[[plot_id_column]]) != "character"){
      stop("`tree_dat[[plot_id_column]]` has to be a `character` vector.")
    }
    if(is.null(tree_dat[[weight_column]])){
      stop("`tree_dat[[weight_column]]` does not exist.")
    }
    if(mode(tree_dat[[weight_column]]) != "numeric"){
      stop("`tree_dat[[weight_column]]` has to be a numeric.")
    }
    if(!all(agg_class_column %in% names(tree_dat))){
      stop("All classes that are specified in `agg_class_column` must be columns of `tree_dat`.")
    }
    if(!all(c(dis_trait_column, recursive = TRUE)  %in% names(tree_dat))){
      stop("All traits that are specified in `dis_trait_column` must be columns of `tree_dat`.")
    }
    
    
    if(is.null(plot_dat[[plot_id_column]])){
      stop("`plot_dat[[plot_id_column]]` does not exist.")
    }
    if(mode(plot_dat[[plot_id_column]]) != "character"){
      stop("`plot_dat[[plot_id_column]]` has to be a `character` vector.")
    }
    core_column <- force(core_column)
    plot_id_values <- check_common_relations(plot_dat, tree_dat, plot_id_column)
    tree_dat <- subset(tree_dat, tree_dat[[plot_id_column]] %in% plot_id_values)
    plot_dat <- subset(plot_dat, plot_dat[[plot_id_column]] %in% plot_id_values)
    plot_dat <- sf::st_set_geometry(plot_dat, core_column)
    apa_map_column <-
      get_voronoi_list(
        core = plot_dat,
        point = tree_dat,
        plot_id_column = plot_id_column,
        buffer_column = buffer_column,
        weight_column = weight_column,
        class_columns = tree_id_column,
        res = res)
    plot_dat$apa_map <- apa_map_column
    
    apa_list <- new_apa_list(c(plot_dat = list(plot_dat), tree_dat = list(tree_dat)),
                             apa_config = list(
                               plot_id_column = plot_id_column,
                               tree_id_column = tree_id_column, 
                               core_column = core_column,
                               weight_column = weight_column,
                               buffer_column = buffer_column,
                               randomized = FALSE,
                               apa_polygon = FALSE,
                               res = res))
    apa_list <- apa_add_agg_class(apa_list, agg_class_column = agg_class_column, apa_polygon = FALSE)
    
    if(tree_id_column %in% subplot_id_column){
      stop("id columns of `tree_dat` or `subplot_dat` may not have the same name.")
    }
    apa_list <- apa_add_subplot_dat(apa_list, subplot_dat,
                                    subplot_id_column = subplot_id_column,
                                    radius = radius,
                                    apa_polygon = FALSE)
    
    if(apa_polygon){
      message("\nStarting polygonization:")
      apa_list <- apa_add_polygon(apa_list)
    }
    
    if(!is.na(apa_properties)[1]){
      prop_order <- c("border_tree", "apa_size", "ndiv", "pdiv")
      apa_properties_proc <- intersect(prop_order, apa_properties)
      missing_prop <- setdiff(apa_properties, apa_properties_proc)
      if(length(missing_prop) > 0){
        stop("`", paste(missing_prop, collapse = "`, `", "` is/are no valid `apa_properteis`. Only `border_tree`, `apa_size` or `ndiv` are accepted (or NA for no calculation of properties)."))
      }
      if("border_tree" %in% apa_properties_proc){
        apa_list <- apa_border_tree(apa_list)
      }
      if("apa_size" %in% apa_properties_proc){
        apa_list <- apa_size(apa_list, edge_correction = edge_correction)
      }
      if(any(c("pdiv", "ndiv") %in% apa_properties_proc)){
        apa_list <- apa_ndiv(apa_list, edge_correction = edge_correction, 
                             dis_trait_column = dis_trait_column,
                             dis_method = dis_method,
                             dis_transform = dis_transform,
                             scope = scope,
                             pdiv = ifelse("pdiv" %in% apa_properties_proc, TRUE, FALSE))
      }
    }
    apa_list
    
  }

#' @keywords internal 
new_apa_list <-
  function(apa_list,
           add_agg_to_apa_list = NULL,
           add_subplot_to_apa_list = NULL,
           apa_config = NULL){
    apa_config_old <- attr(apa_list, "apa_config")
    if(!is.null(add_agg_to_apa_list)){
      apa_list <- c(apa_list, add_agg_to_apa_list)
    }
    if(!is.null(add_subplot_to_apa_list)){
      apa_list$subplot_dat <- c(apa_list$subplot_dat, add_subplot_to_apa_list)
    }
    apa_config$plot_id_values <- apa_list$plot_dat[[apa_config$plot_id_column]]
    apa_config <- update_apa_config(apa_list, apa_config)
    attr(apa_list, "apa_config") <- apa_config
    
    class(apa_list) <- c("apa_list", "list")
    validate_apa_list(apa_list)
  }

#' @keywords internal 
validate_apa_list <- function(apa_list){
  force(apa_list)
  apa_config <- attr(apa_list, "apa_config")
  
  if(!inherits(apa_list, "list")){
    stop("`apa_list` must be a `list`-object.")
  }
  
  if(!names(apa_list)[1] == "plot_dat"){
    stop("The first two elements in `apa_list` must be names `plot_dat` and `tree_dat`")
  }
  
  apa_map_column <- apa_list$plot_dat$apa_map
  if(is.null(apa_map_column)){
    stop("There is no `apa_map`-column in `apa_map$plot_dat`")
  }
  if(is.null(apa_list$plot_dat[[apa_config$plot_id_column]])){
    stop("There is no `plot_id_column` in `apa_list$plot_dat")
  }
  if(!all.equal(names(apa_map_column), apa_config$plot_id_values)){
    stop("Names of `apa_list$plot_dat$apa_map_column` do not match with `apa_config$plot_id_values`")
  }
  if(!all.equal(apa_list$plot_dat[[apa_config$plot_id_column]], apa_config$plot_id_values)){
    stop("Names of `apa_list$plot_dat$apa_map_column` do not match with `apa_config$plot_id_values`")
  }
  if(any(duplicated(apa_config$plot_id_values))){
    stop("`plot_id_column` has to be a unique identifier of `apa_map$plot_dat`")
  }
  apa_map_layer_names <- lapply(apa_map_column, names)
  if(!all(sapply(apa_map_layer_names, length) == length(apa_map_layer_names[[1]]))){
    stop("All `RasterStack` objects in `apa_list$plot_dat$apa_map` need to have identical layer numbers and names.")
  }
  if(nrow(unique(as.data.frame(do.call(rbind, apa_map_layer_names)))) > 2){
    stop("All `RasterStack` objects in `apa_list$plot_dat$apa_map` need to have identical layer numbers and names.")
  }
  if(nrow(unique(as.data.frame(t(sapply(lapply(apa_map_column, raster::res), round, 5))))) > 1){
    stop("All `RasterStack` objects in `apa_list$plot_dat$apa_map` need to have the same resolution.")
  }
  if(nrow(unique(as.data.frame(sapply(lapply(apa_map_column, raster::res), round, 5)))) > 1){
    stop("All `RasterStack` objects in `apa_list$plot_dat$apa_map` need to have identical resolutions in `x`- and `y`-direction.")
  }
  apa_map_layer_names <- apa_map_layer_names[[1]]
  if(!apa_map_layer_names[1] == apa_config$tree_id_column){
    stop("The first layer in the `RasterStack`-objects in `apa_list$plot_dat$apa_map` has to be identical to `tree_id_column`.")
  }
  if(!utils::tail(apa_map_layer_names, 1) == "critical"){
    stop("The last layer in the `RasterStack`-objects in `apa_list$plot_dat$apa_map` has to have the name `critical`.")
  }
  apa_map_agg_layer_names <- apa_map_layer_names[c(-1, -length(apa_map_layer_names))]
  if(length(apa_map_agg_layer_names) == 0){
    apa_map_agg_layer_names <- NA
  }
  if((is.na(apa_config$agg_class_column[1]) & !is.na(apa_map_agg_layer_names[1]))||
     !all.equal(apa_map_agg_layer_names, apa_config$agg_class_column)){
    stop("Layer names of the `RasterStack`-objects in `apa_list$plot_dat$apa_map` have to contain all variables listed in `agg_class_column`.")
  }
  
  stand_col_names <- lapply((apa_list[!grepl("subplot_dat", names(apa_list))]), names)
  subplot_dat_col_names <- lapply(apa_list$subplot_dat, lapply, names)
  col_names <- c(list(stand = stand_col_names), subplot_dat_col_names)
  
  #check naming structure of `apa_list`
  if(length(subplot_dat_col_names) > 0 && !all.equal(as.vector(sapply(lapply(subplot_dat_col_names, names), `[`, 1)),
                                                     names(apa_config$subplot_id_column))){
    stop("The first elements of the lists in `subplot_dat` have to match with `subplot_id_column`.")
  }
  if(!all(sapply(lapply(col_names, names), `[`, 2) == "tree_dat")){
    stop("The second element of `apa_list` and the the second element of all lists in `apa_list$subplot_dat` (if existing) has to have the name `tree_dat`.")
  }
  if(!is.na(apa_config$agg_class_column)[1]){
    if(!all(sapply(lapply(lapply(col_names, names), `%in%`, x = apa_config$agg_class_column), all))){
      stop("All variables specified in `agg_class_column` have to be elements of `apa_list` and of all lists in `apa_list$subplot_dat` (if existing).")
    }
  }
  
  # check plot_id_column and all id_columns
  if(!all(rapply(col_names, f = `%in%`, classes = "character",
                 apa_config$plot_id_column, how = "unlist",
                 x = apa_config$plot_id_column))){
    stop("All datasets in `apa_list` have to have a `plot_id_column`.")
  }
  apa_list_reordered <- 
    c(list(stand = apa_list[!grepl("subplot_dat", names(apa_list))]),
      apa_list$subplot_dat)
  plot_id_column_mode <- 
    c(lapply_deep(apa_list_reordered, function(x){mode(x[[apa_config$plot_id_column]])}), recursive = TRUE)
  if(!all(plot_id_column_mode == "character")){
    stop("`plot_id_column` has to be a `character` vector in all data-sets of `apa_list`.")
  }
  if(!all(c(lapply(apa_list_reordered, lapply, `[[`, apa_config$plot_id_column), recursive = TRUE) %in% apa_config$plot_id_values)){
    stop("All observations in `plot_id_column`s of all datasets have to be values from `apa_config$plot_id_values`.")
  }
  if(!all(sapply(apa_list_reordered, sapply, inherits, "data.frame"))){
    stop("All datasets in `apa_list` have to be `data-frame`-like objects.")
  }
  if(!all(sapply(sapply(col_names, `[`, "tree_dat"), `%in%`, x = apa_config$tree_id_column))){
    stop("`apa_list$tree_dat` and all `tree_dat` datasets in `apa_list$subplot_dat` ", 
         "(if existing) have to have a `tree_id_column`.")
  }
  tree_id <- apa_list$tree_dat[[apa_config$tree_id_column]]
  if(mode(tree_id) != "character"){
    stop("`apa_list$tree_dat[[tree_id_column]]` has to be a `character` vector.")
  }
  if(any(duplicated(tree_id))){
    stop("`apa_list$tree_dat[[tree_id_column]]` has to be a unique identifier of `apa_list$tree_dat`.")
  }
  for(subplot_i in names(apa_config$subplot_id_column)){
    if(!all(sapply(col_names[[subplot_i]], `%in%`, x = apa_config$subplot_id_column[subplot_i]))){
      stop("All datasets in `subplot_dat` have to have id_columns as specified in `subplot_id_columns.")
    }
    if(any(duplicated(apa_list$subplot_dat[[subplot_i]][[subplot_i]][[apa_config$subplot_id_column[subplot_i]]]))){
      stop("`subplot_id_column` has to point to unique identifiers.")
    }
  }
  
  # check existing `apa_polygon`-column
  if(apa_config$apa_polygon && !all(sapply(apa_list_reordered, sapply, inherits, "sf"))){
    stop("If `apa_polygon == TRUE`, all datasets in `apa_list` have to be `sf`-objects")
  }
  if(apa_config$apa_polygon & !is.na(apa_config$agg_class_column)[1]){
    agg_class_data <- lapply(apa_list_reordered, `[`, -1:-2)
    agg_sf_columns <- 
      lapply_deep(agg_class_data, .what = "sf", attr, which = "sf_column")
    agg_sf_columns <- c(agg_sf_columns, recursive = TRUE)
    if(!all(agg_sf_columns == "apa_polygon")){
      stop("If `apa_polygon == TRUE`, apa_polygon` must be the `sf_column` of all aggregated datasets in `apa_list`.")
    }
  }
  
  # check if apa_properties are existing
  if("apa_size" %in% apa_config$apa_properties){
    if(!all(rapply(lapply(col_names, `[`, -1), classes = "character", `%in%`, x = "apa_size"))){
      stop("If `apa_size` is specified in `apa_properties`, all datasets have to have `apa_size`-columns (except for the `plot_dat` and `subplot_dat` datasets).")
    }
  }
  if(any(grepl("_ndiv", apa_config$apa_properties))){
    ndiv_names <- paste0(names(apa_config$dis_trait_column), "_ndiv")
    for(ndiv_name_i in ndiv_names){
      if(!all(rapply(col_names, classes = "character", `%in%`, x = ndiv_name_i))){
        stop("If `ndiv` is specified in `apa_properties`, all datasets have to have `ndiv`-columns.")
      }
    }
  }
  if(any(grepl("_pdiv", apa_config$apa_properties))){
    pdiv_names <- paste0(names(apa_config$dis_trait_column), "_pdiv")
    for(pdiv_name_i in pdiv_names){
      if(!all(rapply(col_names, classes = "character", `%in%`, x = pdiv_name_i))){
        stop("If `ndiv` is specified in `apa_properties`, all datasets have to have `ndiv`-columns.")
      }
    }
  }
  if("border_tree" %in% apa_config$apa_properties){
    if(!all(rapply(lapply(col_names, `[`, "tree_dat"), classes = "character", `%in%`, x = "border_tree"))){
      stop("If `border_tree` is specified in `apa_properties`, all datasets have to have `border_tree`-columns.")
    }
  }
  if(!inherits(apa_list, "apa_list")){
    stop("Class is not `apa_list`.")
  }
  apa_list
}

#' Subsetting of apa_lists
#'
#' Select one or more APA-maps and related data from an \code{apa_list} object.
#'
#' @param x A /code{apa_list}-object.
#' 
#' @param subset Either a numeric or a character vector specifying the plots that
#'   shall be selected.
#'   
#' @param ... not implemented.
#'
#' @method subset apa_list
#' @export
subset.apa_list <- function(x, subset, ...){
  if(is.null(subset)){
    return(x)
  }
  apa_config <- attr(x, "apa_config")
  if(is.numeric(subset)){
    subset <- apa_config$plot_id_values[subset]
  }
  if(!is.character(subset)){
    stop("Invalid `subset`.")
  }
  if(!all(subset %in% apa_config$plot_id_values)){
    stop("All values in `subset` have to be elements of `apa_config$plot_id_column`.")
  }
  apa_stand <- x[c("plot_dat", "tree_dat", apa_config$agg_class_column)]
  for(data_i in grep("subplot_dat", names(x), value = TRUE, invert = TRUE)){
    select_idx <- x[[data_i]][[apa_config$plot_id_column]] %in% subset
    x[[data_i]] <- x[[data_i]][select_idx ,]
  }
  
  for(subplot_j in names(x$subplot_dat)){
    for(data_i in names(x$subplot_dat[[subplot_j]])){
      select_idx <- x$subplot[[subplot_j]][[data_i]][[apa_config$plot_id_column]] %in% subset
      x$subplot_dat[[subplot_j]][[data_i]] <- 
        x$subplot_dat[[subplot_j]][[data_i]][select_idx ,]
    }
  }
  apa_config_new <- apa_config
  apa_config_new$plot_id_values <- subset
  apa_config_new <- update_apa_config(x, apa_config_new)
  attr(x, "apa_config") <- apa_config_new
  validate_apa_list(x)
}

#' Add aggregation classes to an \code{apa_list}
#'
#' Additional grouping variables are added to an \code{apa_list} (for example
#' tree species information). All APA-properties will be aggregated for these
#' classes as well.
#'
#' @template apa_list_arg
#' 
#' @param agg_class_column A vector of column names in \code{apa_list$tree_dat}
#'   that specify variables that are used as additional grouping variable.
#' 
#' @param apa_polygon logical, specifies if \code{POLYGON}s of APA-patches be
#'   should be added.
#' 
#' 
#' @export
apa_add_agg_class <- 
  function(apa_list, agg_class_column = NULL,
           apa_polygon = attr(apa_list, "apa_config")$apa_polygon){
    force(apa_polygon)
    apa_list <- apa_drop_polygon(apa_list)
    if(is.null(agg_class_column) || is.na(agg_class_column[1])){
      return(apa_list)
    }
    apa_config_old <- attr(apa_list, "apa_config")
    if(!is.na(apa_config_old$apa_properties[1])){
      if(any(c("apa_size", "ndiv") %in% apa_config_old$properties)){
        stop("Addition of aggregation classes is only possible before calculating `apa_size` or `ndiv`.")
      }
    }
    if(!is.na(apa_config_old$subplot_id_column[1])){
      stop("Addition of aggregation classes is only possible before adding subplots to `apa_list`.")
    }
    if(!all(agg_class_column %in% names(apa_list$tree_dat))){
      stop("`", paste(agg_class_column, collapse = "`, `"),
           "` is/are no column/s of `tree_dat`.")
    }
    check_agg_class_duplicate <- intersect(agg_class_column, 
                                           apa_config_old$agg_class_column)
    if(length(check_agg_class_duplicate) > 0){
      stop("`", paste(check_agg_class_duplicate, collapse = "`, `"),
           "` is/are already datasets in `apa_list`.")
      
    }
    
    apa_map_column <- apa_list$plot_dat$apa_map
    agg_class_add <- setdiff(agg_class_column,apa_config_old$agg_class_column)
    if(length(agg_class_add) > 0){
      reclass_df <- sf::st_set_geometry(apa_list$tree_dat, NULL)
      reclass_df <- 
        reclass_df[, c(apa_config_old$tree_id_column, agg_class_column), drop = FALSE]
      apa_map_column <- lapply(apa_map_column, rst_reclassify_chr, reclass_df)
      apa_list$plot_dat$apa_map <- apa_map_column
    }
    agg_dat <- names_to_list(agg_class_column)
    message("\nAdding aggregation classes:")
    for(class_i in agg_class_column){
      message(class_i)
      keep_columns <- c(apa_config_old$plot_id_column, class_i)
      agg_dat[[class_i]] <- 
        sf::st_drop_geometry(apa_list$tree_dat)[, keep_columns, drop = FALSE]
      
      agg_dat[[class_i]] <- unique(agg_dat[[class_i]])
    }
    apa_config_new <- apa_config_old
    apa_config_new$agg_class_column <- 
      as.vector(stats::na.omit(c(apa_config_old$agg_class_column,
                                 agg_class_column)))
    apa_list <- new_apa_list(apa_list, 
                             add_agg_to_apa_list = agg_dat,
                             apa_config = apa_config_new)
    if(apa_polygon){
      apa_list <- apa_add_polygon(apa_list)
    }
    apa_list
  }

#' Add subplots to datasets
#'
#' If there are subplots within the plots for which APA-characteristics should be
#' calculated, they may be added to the \code{apa_list}.
#'
#' @template apa_list_arg
#'
#' @param subplot_dat A named list with datasets (\code{sf}-\code{data.frame}s)
#'   about the subplots . Each \code{data.frame} in \code{subplot_dat} has to
#'   have a unique name, a column with the \code{plot_id}, an own id-column
#'   (specified via the \code{subplot_id_column} argument) and a geometry
#'   column. The geometry column either has to contain \code{POLYGON}-data or
#'   \code{POINT}-data (in which case APA-properties in circular neighborhoods
#'   around the points will be calculated.)
#'
#' @param subplot_id_column A named vector specifying the id-columns of the
#'   subplot datasets.
#'
#' @param radius If a geometry column in \code{subplot_dat} contains
#'   \code{POINT}-data, \code{radius} specifies the radius of the neighborhood
#'   analysis.
#'
#' @param apa_polygon Should polygons of the the APA-patches added to the
#'   dataset?
#'
#' @export
#' 
apa_add_subplot_dat <- 
  function(apa_list, 
           subplot_dat, 
           subplot_id_column,
           radius = NULL,
           apa_polygon = attr(apa_list, "apa_config")$apa_polygon){
    if(length(subplot_dat) == 0){
      return(apa_list)
    }
    force(apa_polygon)
    apa_list <- apa_drop_polygon(apa_list)
    apa_config_old <- attr(apa_list, "apa_config")
    plot_id_column <- apa_config_old$plot_id_column
    if(!is.na(apa_config_old$apa_properties[1])){
      stop("Addition of subplots is only possible before calculating APA-properties. 
         Use `apa_drop_properties()` first.")
    }
    if(is.null(radius) || is.na(radius)){
      radius <- apa_config_old$radius
    }
    
    check_subplot_exists <-
      is.logical(all.equal(sort(names(subplot_dat)), 
                           sort(names(subplot_id_column))))
    if(!check_subplot_exists){
      stop("Names of datasets in `subplot_dat` and names of `subplot_id_column`",
           "have to be identical.")
    }
    subplot_columns <- lapply(subplot_dat, names)
    check_plot_id_column<-
      c(lapply(subplot_columns, FUN = `%in%`, x = plot_id_column), recursive = TRUE)
    if(!all(check_plot_id_column)){
      stop("`plot_id_column` (" , plot_id_column, ") has to be present in all subplot datasets.\n")
    }
    
    check_subplot_sf <- sapply(subplot_dat, inherits, "sf")
    if(!all(check_subplot_sf)){
      stop("All elements in `subplot_dat` have to be `sf`-objects.")
    }
    
    check_subplot_point <- 
      any(sapply(lapply(subplot_dat, sf::st_geometry), inherits, "sfc_POINT"))
    if(check_subplot_point & is.na(radius)){
      stop("`radius` has to be specified if a `subplot`` dataset contains point data.")
    }
    if(!is.na(apa_config_old$radius) && radius != apa_config_old$radius){
      stop("One `apa_list` can't have a different `radius` for several point datasets in `subplot_dat` (not yet implemented).")
    }
    
    if(!is.null(apa_config_old$subplot_id_column)){
      check_subplot_id_duplicates <- 
        intersect(names(apa_config_old$subplot_id_column),
                  names(subplot_id_column))
      if(length(check_subplot_id_duplicates) > 0){
        stop("`", paste(check_subplot_id_duplicates, collapse = "`, `"),
             "` is/are already subplot-dataset in `apa_list`.")
      }
    }
    output_list <- names_to_list(names(subplot_id_column))
    message("\nAdding subplot datasets:")
    for(subplot_i in names(subplot_id_column)){
      message("`", subplot_i, "` - setting up datasets")
      if(!subplot_id_column[subplot_i] %in% names(subplot_dat[[subplot_i]])){
        stop("`", subplot_id_column[subplot_i], "' is not a column of `subplot_dat$",
             subplot_i, "`,")
      }
      subplot_dat[[subplot_i]] <- 
        subset(subplot_dat[[subplot_i]], subplot_dat[[subplot_i]][[plot_id_column]] %in% apa_config_old$plot_id_values)
      output_list[[subplot_i]] <- list()
      output_list[[subplot_i]][[subplot_i]] <- subplot_dat[[subplot_i]]
      subplot_intersect <- subplot_dat[[subplot_i]][c(plot_id_column, subplot_id_column[subplot_i])]
      if(is.null(subplot_intersect)){
        stop("'", subplot_i, "' not found in 'subplot_dat'.")
      }
      if(inherits(sf::st_geometry(subplot_intersect), "sfc_POINT")){
        if(is.null(radius) || is.na(radius)){
          stop("Please specify neighborhood radius for point data in `subplot_dat` via `apa_config = list(radius == xxx)`.")
        }
      }
      sf::st_agr(subplot_intersect) <- "constant"
      output_list[[subplot_i]]$tree_dat <- 
        class_area(apa_list$plot_dat$apa_map, subplot = subplot_intersect,
                   plot_id_column = plot_id_column,
                   layer = stats::na.omit(c(apa_config_old$tree_id_column,
                                            apa_config_old$agg_class_column)),
                   radius = radius)
      match_idx <- match_by(output_list[[subplot_i]]$tree_dat, apa_list$tree_dat, by = apa_config_old$tree_id_column)
      subplot_tree_geo <- sf::st_geometry(apa_list$tree_dat)[match_idx]
      output_list[[subplot_i]]$tree_dat[[attr(apa_list$tree_dat, "sf_column")]] <- 
        subplot_tree_geo
      output_list[[subplot_i]]$tree_dat <- sf::st_as_sf(output_list[[subplot_i]]$tree_dat)
      output_list[[subplot_i]]$tree_dat$area <- NULL
      if(!is.na(apa_config_old$agg_class_column)[1]){
        for(class_i in apa_config_old$agg_class_column){
          agg_class_data <- apa_list[[class_i]]
          if(inherits(agg_class_data, "sf")){
            agg_class_data <- sf::st_drop_geometry(agg_class_data)
          }
          agg_class_list <- split(agg_class_data, agg_class_data[[plot_id_column]])
          agg_class_list <- lapply(agg_class_list, `[`, class_i)
          subplot_intersect_no_sf <- sf::st_drop_geometry(subplot_intersect)
          subplot_intersect_list <- split(subplot_intersect_no_sf, subplot_intersect_no_sf[[plot_id_column]])
          subplot_intersect_list <- lapply(subplot_intersect_list, lapply, unique)
          common_relations <- check_common_relations(subplot_intersect_list, agg_class_list)
          joined_list <- mapply(FUN = c, subplot_intersect_list[common_relations],
                                agg_class_list[common_relations])
          full_combinations <- apply(joined_list, 2, expand.grid, stringsAsFactors = FALSE)
          full_combinations <- do.call("rbind", full_combinations)
          full_combinations <- full_combinations[do.call(order, full_combinations), ]
          
          output_list[[subplot_i]][[class_i]] <- full_combinations
          
        }
      }
    }
    message("")
    apa_config_new <- apa_config_old
    apa_config_new$radius = radius
    apa_config_new$subplot_id_column <- c(apa_config_old$subplot_id_column, subplot_id_column)
    apa_config_new$subplot_id_column <- apa_config_new$subplot_id_column[!is.na(apa_config_new$subplot_id_column)]
    
    apa_list <- 
      new_apa_list(apa_list, add_subplot_to_apa_list = output_list, 
                   apa_config = apa_config_new)
    if(apa_polygon){
      apa_list <- apa_add_polygon(apa_list)
    }
    apa_list
  }

#' Add \code{POLYGON}-columns to an \code{apa_list}
#'
#' @template apa_list_arg
#'
#' @export
apa_add_polygon <- function(apa_list){
  apa_config <- attr(apa_list, "apa_config")
  apa_map_column <- apa_list$plot_dat$apa_map
  message("`tree_dat` - polygonizing APA-patches:")
  tree_apa_polygon <-
    polygonize_class(apa_map_column, layer = apa_config$tree_id_column,
                     plot_id_column = apa_config$plot_id_column)
  match_tree_id <- match(apa_list$tree_dat[[apa_config$tree_id_column]],
                         tree_apa_polygon[[apa_config$tree_id_column]])
  apa_polygon_column <- 
    sf::st_sfc(rep(list(sf::st_polygon()), nrow(apa_list$tree_dat)),
               crs = sf::st_crs(apa_list$tree_dat))
  apa_polygon_column[!is.na(match_tree_id)] <- tree_apa_polygon[[3]][stats::na.omit(match_tree_id)]
  apa_list$tree_dat$apa_polygon <- apa_polygon_column[1:length(apa_polygon_column)]
  sf::st_agr(apa_list$tree_dat) <- "constant"
  
  if(!is.na(apa_config$agg_class_column)[1]){
    for(class_i in apa_config$agg_class_column){
      message("`", class_i, "` - polygonizing APA-patches:  ", appendLF = FALSE)
      agg_polygon_i <- 
        polygonize_class(rst = apa_map_column, layer = class_i,
                         plot_id_column = apa_config$plot_id_column,
                         geometry_name = "apa_polygon")
      match_idx <- match_by(apa_list[[class_i]], agg_polygon_i,
                            by = c(apa_config$plot_id_column, class_i))
      apa_list[[class_i]]$apa_polygon <- agg_polygon_i[match_idx, ]$apa_polygon
      apa_list[[class_i]] <- sf::st_sf(apa_list[[class_i]])
      apa_list[[class_i]] <- suppressWarnings(sf::st_set_crs(apa_list[[class_i]], sf::st_crs(apa_list$plot_dat)))
      apa_list[[class_i]] <- suppressWarnings(sf::st_set_agr(apa_list[[class_i]], "constant"))
    }
  }
  
  if(!is.na(apa_config$subplot_id_column[1])){
    for(subplot_i in names(apa_config$subplot_id_column)){
      id_columns <- 
        c(apa_config$plot_id_column, apa_config$subplot_id_column[subplot_i])
      subplot_intersect <- apa_list$subplot_dat[[subplot_i]][[subplot_i]][id_columns]
      if(is.null(subplot_intersect)){
        stop("'", subplot_i, "' not found in 'subplot_dat'.")
      }
      if(inherits(sf::st_geometry(subplot_intersect), "sfc_POINT")){
        if(is.null(apa_config$radius) || is.na(apa_config$radius)){
          stop("Please specify neighborhood radius for point data in `subplot_dat` via `apa_config = list(radius == xxx)`.")
        }
        subplot_intersect <- sf::st_buffer(subplot_intersect, apa_config$radius)
      }
      sf::st_agr(subplot_intersect) <- "constant"
      tree_intersect <- apa_list$tree_dat
      sf::st_geometry(tree_intersect) <- "apa_polygon"
      tree_intersect <- tree_intersect[, c(apa_config$plot_id_column, apa_config$tree_id_column)]
      intersection_apa_polygones <- 
        sf::st_intersection(y = subplot_intersect[, apa_config$subplot_id_column[subplot_i]],
                            tree_intersect)
      join_columns <- c(id_columns, apa_config$tree_id_column)
      match_idx <- match_by(apa_list$subplot_dat[[subplot_i]]$tree_dat,
                            intersection_apa_polygones,
                            by = join_columns)
      apa_list$subplot_dat[[subplot_i]]$tree_dat$apa_polygon <- 
        intersection_apa_polygones$apa_polygon[match_idx]
      apa_list$subplot_dat[[subplot_i]]$tree_dat <-
        sf::st_sf(apa_list$subplot_dat[[subplot_i]]$tree_dat)
      
      if(!is.na(apa_config$agg_class_column)[1]){
        for(class_i in apa_config$agg_class_column){
          intersection_apa_polygones <- 
            sf::st_intersection(apa_list[[class_i]],
                                subplot_intersect[, apa_config$subplot_id_column[subplot_i]])
          join_columns <- c(id_columns, class_i)
          match_idx <- match_by(apa_list$subplot_dat[[subplot_i]][[class_i]],
                                intersection_apa_polygones,
                                by = join_columns)
          apa_list$subplot_dat[[subplot_i]][[class_i]]$apa_polygon <- 
            intersection_apa_polygones$apa_polygon[match_idx]
          apa_list$subplot_dat[[subplot_i]][[class_i]] <-
            sf::st_sf(apa_list$subplot_dat[[subplot_i]][[class_i]])
        }
      }
    }
  }
  apa_config_new <- update_apa_config(apa_list, apa_config = list(apa_polygon = TRUE))
  attr(apa_list, "apa_config") <- apa_config_new
  validate_apa_list(apa_list)
}  

#' @keywords internal 
update_apa_config <- function(apa_list, apa_config){
  apa_config_min <- list(plot_id_column = NA,
                         plot_id_values = NA,
                         tree_id_column = NA,
                         core_column = NA,
                         buffer_column = NA, 
                         weight_column = NA,
                         res = NA,
                         apa_polygon = NA,
                         agg_class_column = NA,
                         subplot_id_column = NA,
                         radius = NA,
                         apa_properties = NA,
                         edge_correction = NA,
                         dis_trait_column = NA,
                         dis_method = NA,
                         dis_transform = NA,
                         scope = NA,
                         randomized = NA)
  apa_config_old <- attr(apa_list, "apa_config")
  apa_config_old <- 
    apa_config_old[!names(apa_config_old) %in% c("names", "class")]
  if(is.null(apa_config)){
    apa_config <- apa_config_old
  }
  
  apa_config_remain <- 
    apa_config_old[setdiff(names(apa_config_old), names(apa_config))]
  apa_config[names(apa_config_remain)] <- apa_config_remain
  
  apa_config_remain <- 
    apa_config_min[setdiff(names(apa_config_min), names(apa_config))]
  apa_config[names(apa_config_remain)] <- apa_config_remain
  apa_config <- apa_config[names(apa_config_min)]
  apa_config[sapply(apa_config, is.null)] <- NA
  apa_config <- apa_config[names(apa_config_min)]
  return(apa_config)
}

#' @keywords internal
update_apa <- function(apa_list, apa_config_new = NULL, force = NULL, no_polygon_update = FALSE){
  if(is.null(apa_config_new) & is.null(force)){
    return(NULL)
  }
  update_support <- 
    c("plot_id_column", "core_column", "buffer_column", "weight_column",
      "apa_properties", "radius", "edge_correction", "dis_trait_column",
      "dis_method", "dis_transform", "scope", "apa_polygon")
  if(any(!names(apa_config_new) %in% update_support)){
    stop("Currently only the following `apa_config` elements can be updated: `",
         paste(update_support, collapse = "`, `"), "`.")
  }
  apa_config_old <- attr(apa_list, "apa_config")
  force_map_update <- FALSE
  force_polygon_update <- FALSE
  force_agg_class_update <- FALSE
  force_subplot_update <- FALSE
  force_properties_update <- FALSE
  if(any(names(apa_config_new) %in% 
         c("plot_id_column", "core_column", "buffer_column", "weight_column"))){
    force <- unique(c(force, "map"))
  }
  if(any(c("all", "map") %in% force)){
    force_map_update <- TRUE
    force <- unique(c(force, "agg_class"))
  }
  if(any(c("agg_class") %in% force)){
    force_agg_class_update <- TRUE
    force <- unique(c(force, "subplot"))
  }
  if("radius" %in% names(apa_config_new)){
    force <- unique(c(force, "subplot"))
  }
  if("subplot" %in% force){
    force <- unique(c(force, "properties"))
    if(apa_config_old$apa_polygon){
      force <- unique(c(force, "polygon"))
    }
    force_subplot_update <- TRUE
  }
  
  if(any(c("apa_properties", "edge_correction", "dis_trait_column",
           "dis_method", "dis_transform", "scope") %in% names(apa_config_new))){
    force <- unique(c(force, "properties"))
  }
  if("properties" %in% force){
    force_properties_update <- TRUE
  }
  
  if(!is.null(apa_config_new$apa_polygon)){
    if(apa_config_new$apa_polygon){
      force <- unique(c(force, "polygon"))}
    else{
      no_polygon_update <- TRUE
    }
  } 
  if("polygon" %in% force & !no_polygon_update){
    force_polygon_update <- TRUE
  }
  
  apa_list_updated <- apa_list
  apa_config_updated <- update_apa_config(apa_list, apa_config_new)
  
  if(force_map_update){
    apa_list_updated <- apa_drop_polygon(apa_list_updated)
    apa_list_updated <- apa_drop_subplot(apa_list_updated)
    apa_list_updated <- apa_drop_agg_class(apa_list_updated)
    apa_list_updated <- apa_drop_properties(apa_list_updated)
    plot_dat <- apa_list_updated$plot_dat
    tree_dat <- apa_list_updated$tree_dat
    plot_id_column <- apa_config_updated$plot_id_column
    plot_id_values <- check_common_relations(plot_dat, tree_dat, 
                                             apa_config_updated$plot_id_column)
    sf::st_geometry(plot_dat) <- apa_config_updated$core_column
    
    tree_dat <- subset(tree_dat, tree_dat[[plot_id_column]] %in% plot_id_values)
    plot_dat <- subset(plot_dat, plot_dat[[plot_id_column]] %in% plot_id_values)
    apa_map_column <-
      get_voronoi_list(
        core = plot_dat,
        point = tree_dat,
        plot_id_column = plot_id_column,
        buffer_column = apa_config_updated$buffer_column,
        weight_column = apa_config_updated$weight_column,
        class_columns = apa_config_updated$tree_id_column,
        res = apa_config_updated$res)
    plot_dat$apa_map <- apa_map_column
    
    apa_list_updated <- new_apa_list(c(plot_dat = list(plot_dat), tree_dat = list(tree_dat)),
                                     apa_config = list(
                                       plot_id_column = plot_id_column,
                                       tree_id_column = apa_config_updated$tree_id_column, 
                                       core_column = apa_config_updated$core_column,
                                       weight_column = apa_config_updated$weight_column,
                                       buffer_column = apa_config_updated$buffer_column,
                                       randomized = FALSE,
                                       apa_polygon = FALSE,
                                       res = apa_config_updated$res))
  }
  if(force_agg_class_update){
    apa_list_updated <- apa_drop_polygon(apa_list_updated)
    apa_list_updated <- apa_drop_subplot(apa_list_updated)
    apa_list_updated <- apa_drop_agg_class(apa_list_updated)
    apa_list_updated <- apa_drop_properties(apa_list_updated)
    message("Initializing aggregation class datasets.")
    apa_list_updated <- apa_add_agg_class(apa_list_updated,
                                          agg_class_column = apa_config_updated$agg_class_column,
                                          apa_polygon = FALSE)
  }
  if(force_subplot_update){
    apa_list_updated <- apa_drop_polygon(apa_list_updated)
    apa_list_updated <- apa_drop_subplot(apa_list_updated)
    apa_list_updated <- apa_drop_properties(apa_list_updated)
    
    message("Initializing subplot datasets.")
    subplot_dat_original <- 
      lapply(apa_drop_properties(apa_list)$subplot_dat, `[[`, 1)
    apa_list_updated <-
      apa_add_subplot_dat(apa_list_updated, subplot_dat_original,
                          subplot_id_column = apa_config_updated$subplot_id_column,
                          radius = apa_config_updated$radius,
                          apa_polygon = FALSE)
  }
  if(force_properties_update){
    apa_list_updated <- apa_drop_properties(apa_list_updated)
    apa_properties <- apa_config_updated$apa_properties
    if("seg_ndiv" %in% apa_properties){
      warning("`seg_ndiv` will not be updated.",
              " Use `apa_seg()` to calculate the segregation of ndiv")
    }
    radius <- apa_config_updated$radius
    edge_correction <- apa_config_updated$edge_correction
    dis_trait_column <- apa_config_updated$dis_trait_column
    dis_method <- apa_config_updated$dis_method
    dis_transform <- apa_config_updated$dis_transform
    scope <- apa_config_updated$scope
    if(!is.na(apa_properties)[1]){
      pdiv_props <- grep(".*_pdiv$", apa_properties, value = TRUE)
      ndiv_props <- grep(".*_ndiv$", apa_properties, value = TRUE)
      size_prop <- grep("^apa_size|critical", apa_properties, value = TRUE)
      border_tree_prop <- grep("^border_tree$", apa_properties, value = TRUE)
      apa_properties_proc <- 
        c(border_tree_prop, size_prop, ndiv_props, pdiv_props)
      missing_prop <- setdiff(apa_properties, apa_properties_proc)
      if(length(missing_prop) > 0){
        stop("`", paste(missing_prop, collapse = "`, `", "` is/are no valid `apa_properteis`. Only `border_tree`, `apa_size` or `ndiv` are accepted (or NA for no calculation of properties)."))
      }
      if("border_tree" %in% apa_properties_proc){
        apa_list_updated <- apa_border_tree(apa_list_updated)
      }
      if("apa_size" %in% apa_properties_proc){
        apa_list_updated <- apa_size(apa_list_updated, edge_correction = edge_correction)
      }
      if(any(grepl("ndiv", apa_properties_proc))){
        apa_list_updated <- apa_ndiv(apa_list_updated, edge_correction = edge_correction, 
                                     dis_trait_column = dis_trait_column,
                                     dis_method = dis_method,
                                     dis_transform = dis_transform,
                                     scope = scope,
                                     pdiv = ifelse(any(grepl("pdiv", apa_properties_proc)), TRUE, FALSE))
      }
    }
  }
  if(force_polygon_update == TRUE){
    apa_list_updated <- apa_add_polygon(apa_list_updated)
  }
  apa_list_updated
}

#' @method print apa_list
#' @export
print.apa_list <- function(x, ...){
  apa_config <- attributes(x)$apa_config
  cat("List with APA-maps of", nrow(x$plot_dat), 
      c("plot", "plots")[min(2, nrow(x$plot_dat))], "and", 
      nrow(x$tree_dat), 
      c("tree", "trees")[min(2, nrow(x$tree_dat))])
  cat(".\n")
  plot_id_string <- shorten_string(apa_config$plot_id_values, 5, 10)
  cat("plot id column:      ", apa_config$plot_id_column, 
      " (\"", paste(plot_id_string, collapse = '" "'), "\"", sep = "")
  if(length(apa_config$plot_id_values) > 5){
    cat(") with", length(apa_config$plot_id_values) - 5, "more plots\n")
  }else{
    cat(")\n")
  }
  bbox <- sf::st_bbox(x$plot_dat)
  bbox_cat <- paste0(names(bbox), ": ", round(bbox, 1))
  cat("bbox:               ", bbox_cat[c(1, 3)], "\n",
      "                   ", bbox_cat[c(2, 4)])
  cat("\n")
  crs_cat <- suppressMessages(utils::capture.output(print(sf::st_geometry(x$plot_dat))))
  crs_cat <- grep("CRS", crs_cat, value = TRUE)
  crs_cat <- sub("CRS:", "CRS:     ", crs_cat)
  cat(crs_cat)
  cat("\n")
  cat("resolution:         ", apa_config$res)
  cat("\n")
  if(!is.na(apa_config$agg_class_column[1])){
    cat("aggregation classes:")
    agg_class_string <- shorten_string(apa_config$agg_class_column, 3, 15)
    cat(" \"", paste(agg_class_string, collapse = '" "'), "\"", sep = "")
    if(length(apa_config$agg_class_column) > 3){
      cat(" with", length(apa_config$agg_class_column) - 3, "more",
          c("class", "classes")[min(2, length(apa_config$agg_class_column) - 3)], "\n")
    }
    cat("\n")
  }
  if(!is.na(apa_config$subplot_id_column[1])){
    subplot_string <- shorten_string(names(apa_config$subplot_id_column), 3, 15)
    cat("subplot datasets:   ")
    cat(" \"", paste(subplot_string, collapse = '" "'), "\"", sep = "")
    if(length(apa_config$subplot_id_column) > 3){
      cat(" with", length(apa_config$subplot_id_column) - 3, "more",
          c("dataset", "datasets")[min(2, length(apa_config$subplot_id_column) - 3)])
    }
    cat( "\n")
  }
  if(!is.na(apa_config$apa_properties[1])){
    apa_properties_string <- shorten_string(apa_config$apa_properties, 3, 15)
    cat("apa properties:      ")
    cat(apa_properties_string, sep = ", ")
    if(length(apa_config$apa_properties) > 3){
      cat(" with", length(apa_config$apa_properties) - 3, "more",
          c("apa property", "apa_properties")[min(2, length(apa_config$apa_properties) - 3)])
    }
    cat("\n")
  }
}
