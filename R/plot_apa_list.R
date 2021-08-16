#' plot APA-maps
#'
#' APA-maps that are stored in an \code{APA-list}-object are being plotted.
#'
#' @param x A \code{apa_list}-object.
#'
#' @param subset The subset of APA-maps in \code{apa_list} that will be plotted
#'   (see [subset.apa_list]).
#'
#' @param color_map If \code{numeric} or \code{character}, this argument
#'   specifies which aggregation class should be used for coloring APA-patches
#'   (see [apa_add_agg_class()]). By supplying a \code{data.frame}, a custom color
#'   scheme can be specified. The first column of the \code{data.frame} needs to
#'   match with a categorical column in \code{tree_dat} and the second column
#'   specifies the colors.
#'
#' @param tree_size_column Column name in \code{tree_dat} that is used to adjust
#'   point sizes of tree locations. If this column is a metric units object, it
#'   is scaled to meters to match the units of the axes.
#'
#' @param tree_size_scale If \code{tree_size_column} is a metric units object,
#'   it is multiplied with \code{tree_size_scale}. Otherwise,
#'   \code{tree_size_scale} specifies the size of the largest point in meters.
#'
#' @param pal If no \code{color_map} is provided, this color palette function is
#'   used for coloring. specified here.
#'
#' @param single_plots Should the APA-maps be plotted as individual plots
#'   (\code{TRUE}) or should all APA-maps be plotted in a single diagram
#'   (\code{FALSE}, default)?
#'
#' @param add_legend If \code{TRUE} (default), a legend will be plotted.
#'
#' @param critical A color to shade the critical area close to the plot borders
#'   that may be influenced by an edge effect.
#'
#' @param add If \code{FALSE} (default), the plot will be standalone, otherwise
#'   it will be added to the currently active plot.
#'
#' @param cex Character expansion factor for labels (legend, titles)
#'
#' @param add_subplot If \code{TRUE}, subplots that were added via
#'   [apa_add_subplot_dat] will be added to the plot (default is \code{FALSE}).
#'
#' @param add_plot_id_values If \code{TRUE} (default), plot id's will be added
#'   as title.
#'   
#' @param ... not implemented.
#'
#' @example examples/example_plot_apa_list.R
#'
#' @method plot apa_list
#' @export
#' @md
#' 
plot.apa_list <- 
  function(x, 
           subset = NULL, 
           color_map = 1,
           tree_size_column = attr(x, "apa_config")$weight_column,
           tree_size_scale = 1,
           pal = sf::sf.colors,
           single_plots = FALSE,
           add_legend = TRUE, 
           critical = "#FF000033",
           add = FALSE, 
           cex = graphics::par("cex"),
           add_subplot = FALSE,
           add_plot_id_values = TRUE,
           ...){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))            
    x <- subset(x, subset)
    old_mfrow_mai <- graphics::par()[c("mfrow", "mai")]
    apa_config <- attr(x, "apa_config")
    if(add){
      add_legend <- FALSE
      if(length(apa_config$plot_id_values) > 1){
        stop("`add = TRUE` is only possible for `apa_lists`` that contain only a single apa_map. Use `subset` to choose one.")
      }
    }
    if(!apa_config$apa_polygon){
      x <- apa_add_polygon(x)
      apa_config <- attr(x, "apa_config")
    }
    tree_dat_proc <- x$tree_dat
    plot_dat_proc <- x$plot_dat
    if(!is.na(apa_config$subplot_id_column[1])){
      subplot_dat_proc <- lapply(x$subplot_dat, `[[`, 1)
    }
    if(is.numeric(color_map)){
      if(is.na(apa_config$agg_class_column[1])){
        if(is.null(formals(pal)$categorical)){
          tree_dat_proc$.color <- pal(n = nrow(x$tree_dat))
        }else{
          tree_dat_proc$.color <- pal(n = nrow(x$tree_dat), categorical = TRUE)
        }
        color_map <- NULL
      }else{
        color_map <- apa_config$agg_class_column[color_map]
      }
    }
    
    if(is.character(color_map)){
      color_map_dat <- x[[color_map]]
      if(inherits(color_map_dat, "sf")){
        color_map_dat <- sf::st_drop_geometry(color_map_dat)
      }
      if("apa_size" %in% apa_config$apa_properties){
        agg_formula <- stats::as.formula(paste0("apa_size ~ ", color_map))
        color_map_dat <- stats::aggregate(agg_formula, color_map_dat, sum)
        color_map_dat <- 
          color_map_dat[sort(color_map_dat$apa_size, decreasing = TRUE, index.return = TRUE)$ix, ]
        color_map_dat$apa_size <- NULL
      }else{
        color_map_dat <- unique(color_map_dat[color_map])
      }
      if(is.null(formals(pal)$categorical)){
        color_map_dat$.color <- pal(n = nrow(color_map_dat))
      }else{
        color_map_dat$.color <- pal(n = nrow(color_map_dat), categorical = TRUE)
      }
      color_map <- color_map_dat
    }
    
    if(!is.null(color_map)){
      names(color_map)[2] <- ".color"
      match_idx <- match_by(tree_dat_proc, color_map, by = names(color_map)[1])
      tree_dat_proc$.color <- color_map[match_idx, ][[2]]
    }
    
    apa_config <- attr(x, "apa_config")
    plot_id_values <- apa_config$plot_id_values
    if(single_plots){
      for(relation_i in apa_config$plot_id_values){
        arguments_i <- 
          replace(as.list(match.call())[-1],
                  c("subset", "single_plots"), 
                  c(list(subset = relation_i), list(single_plots = FALSE)))
        do.call(plot.apa_list, arguments_i)
      }
      return(invisible(NULL))
    }
    
    max_n <- length(plot_id_values)
    if(!is.null(color_map) & add_legend){
      max_n <- max_n + 1
    }
    row_height_in <- graphics::par("ps") * cex / 72
    if(!add){
      asp <- graphics::par("pin")[1]/graphics::par("pin")[2]
      if(asp < 0){
        stop("Plot size too small, please increase plotting area.")
      }
      ncol <- 1:max_n
      nrow <- floor(ncol/ifelse(add_plot_id_values, 1.1, 1)/asp)
      nrow <- ifelse(nrow == 0, 1, nrow)
      total <- ncol*nrow
      ncol <- which(total >= max_n)[1]
      nrow <- ceiling(max_n/ncol)
      total <- ncol * nrow
      graphics::layout(matrix(1:total, ncol = ncol, byrow = TRUE))
    }
    graphics::par(cex = cex)
    if(!is.null(color_map) & add_legend){
      class_present <- unique(tree_dat_proc[[names(color_map)[1]]])
      color_map_sub <- subset(color_map, color_map[[1]] %in% class_present)
      if(!is.null(x[[names(color_map)[1]]]$apa_size) & length(apa_config$plot_id_values) == 1){
        size_idx <- sort(x[[names(color_map)[1]]]$apa_size, index.return = TRUE, decreasing = TRUE)$ix
        size_order <- x[[names(color_map)[1]]][[names(color_map)[[1]]]][size_idx]
        color_idx <- match(size_order, color_map_sub[[1]])
        color_map_sub <- 
          color_map_sub[color_idx, ]
      }
      graphics::par(mai = c(0, 0, 0, 0))
      legend_height <- graphics::par("pin")[2]
      legend_no_of_rows <- floor(legend_height / row_height_in)
      entry_no <- nrow(color_map_sub)
      if(entry_no > legend_no_of_rows){
        warning("Legend size is too small to show all values. Increase plot size or decrease text size.")
      }
      entry_no <- min(legend_no_of_rows, entry_no)
      row_height_rel <- 1 / entry_no
      graphics::plot.new()
      graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      entry_y <- seq(1 - .5 * row_height_rel, 0 + .5 * row_height_rel, by = -row_height_rel)
      graphics::symbols(x = rep(.1, entry_no), y = entry_y, 
                        circles = rep(1 / entry_no, entry_no),
                        bg = color_map_sub$.color[1:entry_no], add = TRUE,
                        inches = legend_height / legend_no_of_rows / 2)
      if(length(color_map_sub) > 2){
        label <- color_map_sub[[3]]
      }else{
        label <- color_map_sub[[1]]
      }
      graphics::text(x = rep(.2, entry_no), y = entry_y,
                     label = label[1:entry_no], adj = 0)
      max_n <- max_n - 1
    }
    graphics::par(mai = c(0, 0, row_height_in * 1.2 * add_plot_id_values, 0))
    for(relation_i in plot_id_values[1:max_n]){
      tree_dat_i <- 
        subset(tree_dat_proc,
               tree_dat_proc[[apa_config$plot_id_column]] == relation_i)
      plot_dat_i <- 
        subset(plot_dat_proc, 
               plot_dat_proc[[apa_config$plot_id_column]] %in% relation_i)
      plot_dat_i <- 
        subset(plot_dat_proc, 
               plot_dat_proc[[apa_config$plot_id_column]] %in% relation_i)
      bb <- sf::st_union(sf::st_as_sfc(sf::st_bbox(tree_dat_i)),
                         sf::st_as_sfc(sf::st_bbox(tree_dat_i$apa_polygon)))
      bb <- sf::st_union(bb, plot_dat_i[[apa_config$core_column]])
      bb <- sf::st_union(bb, plot_dat_i[[apa_config$buffer_column]])
      plot(bb, border = NA, add = add)
      if(add_plot_id_values){
        graphics::mtext(relation_i, side = 3, line = 0, cex = graphics::par("cex"), padj = .3)
      }
      plot(tree_dat_i$apa_polygon, col = tree_dat_i$.color, add = TRUE, border = NA)
      if(!is.null(critical)){
        apa_map <- x$plot_dat$apa_map[[relation_i]]
        critical_poly <- polygonize_class(apa_map, layer = "critical")
        critical_poly <- sf::st_geometry(critical_poly)[critical_poly$critical == 1]
        plot(critical_poly, col = critical, add = TRUE, border = NA)
      }
      plot(tree_dat_i$apa_polygon[], add = TRUE)
      
      if(!is.null(tree_size_column)){
        if(inherits(tree_dat_i[[tree_size_column]], "units")){
          radii <- tryCatch(units::set_units(tree_dat_i[[tree_size_column]], "m"), error = identity)
          if(inherits(radii, "error")){
            radii <- units::drop_units(tree_dat_i[[tree_size_column]])
          }
        }else{
          radii <- tree_dat_i[[tree_size_column]]
        }
        if(is.null(tree_size_scale)){
          inches <- FALSE
        }else{
          if(inherits(radii, "units")){
            inches <- FALSE
            radii <- radii * tree_size_scale
          }else{
            inches <- graphics::par("pin")[1]/diff(graphics::par("usr")[1:2]) * tree_size_scale
          }
        }
        
        graphics::symbols(sf::st_coordinates(tree_dat_i),
                          circles = radii,
                          inches = inches, add = TRUE, bg = "white")
      }
      plot(plot_dat_i[[apa_config$buffer_column]], add = TRUE, lwd = 2)
      plot(plot_dat_i[[apa_config$core_column]], add = TRUE, lwd = 2)
      if(!is.na(apa_config$subplot_id_column[1]) & add_subplot){
        subplot_dat_i <- 
          lapply(subplot_dat_proc, function(x){subset(x, x[[apa_config$plot_id_column]] %in% relation_i)})
        for(subplot_j in names(apa_config$subplot_id_column)){
          subplot_ij_geo <- sf::st_geometry(subplot_dat_i[[subplot_j]])
          if(inherits(subplot_ij_geo, "sfc_POINT")){
            plot(subplot_ij_geo, add = TRUE, pch = 10)
            subplot_ij_geo <- sf::st_buffer(subplot_ij_geo, apa_config$radius)
          }
          plot(subplot_ij_geo, add = TRUE, lwd = 2)
        }
      }
    }
    if(!add){
      graphics::par(old_mfrow_mai)
    }
    invisible(NULL)
  }
