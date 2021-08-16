#' @keywords internal
get_voronoi_list <- 
  function(core,
           point,
           plot_id_column,
           buffer_column = NULL,
           weight_column =  NULL,
           class_columns = NULL,
           res = 1){
    
    # check input classes
    if(!inherits(core, "sf")){
      stop("'core' is not an 'sf'-object.")
    }
    if(!inherits(point, "sf") || !inherits(sf::st_geometry(point), "sfc_POINT")){
      stop("'point' is not an 'sf'-object with a POINT geometry column.")
    }
    
    # assign core geometry to buffer, if buffer empty (core and buffer identical)
    if(is.null(buffer_column)){
      core$.buffer <- sf::st_geometry(core)
      buffer_column <- ".buffer"
    }
    
    # check if all sepcified columns exist
    if(is.null(core[[plot_id_column]])){
      stop("'", plot_id_column, "' is not a column of 'core'.")
    }
    if(is.null(core[[buffer_column]])){
      stop("'", buffer_column, "' is not a column of 'core'.")
    }
    if(is.null(point[[plot_id_column]])){
      stop("'", plot_id_column, "' is not a column of 'point'.")
    }
    if(!all(class_columns %in% names(point))){
      stop("'", paste(setdiff(class_columns, names(point)), collapse = "', '"),
           "' is not a column of 'point'.")
    }
    
    # check uniqueness of 'plot_id_column' in 'core'
    if(any(table(core[[plot_id_column]]) > 2)){
      stop("'plot_id_column' points towards a column in 'core' that is not unique.")
    }
    
    # check empty features
    if(any(sf::st_is_empty(sf::st_geometry(core)))){
      stop("'core' contains an empty feature")
    }
    if(any(sf::st_is_empty(core[[buffer_column]]))){
      stop("'buffer' contains an empty feature")
    }
    if(any(sf::st_is_empty(point))){
      stop("'point' contains empty features.")
    }
    plot_id_values <- check_common_relations(core, point, plot_id_column)
    
    # Generating single APA-maps and store them in a list
    output_list <- names_to_list(plot_id_values)
    message("Generating APA-maps:")
    for(relation_i in plot_id_values){
      message(relation_i, "  ", appendLF = FALSE)
      core_subset <- subset(core, core[[plot_id_column]] == relation_i)
      single_core <- sf::st_geometry(core_subset)
      single_buffer <- core_subset[[buffer_column]]
      point_subset <- subset(point, point[[plot_id_column]] == relation_i)
      output_list[[relation_i]] <- 
        get_single_voronoi(
          core =  single_core,
          point = point_subset,
          buffer = single_buffer,
          weight_column = weight_column,
          class_columns = class_columns,
          res = res,
          warn = FALSE)
    }
    message("")
    output_list
  }
