#' Calculate the segregation index of NDiv
#'
#' The segregation index of NDiv assesses how much NDiv deviates from its random
#' expectation using a randomization approach. Values close to zero indicate a
#' random distribution of the traits that were used to calculate NDiv, values
#' larger than zero indicate a regular distribution (stem-wise mixing) and
#' values smaller than zero indicate a clustered configuration (patch-wise
#' mixing).
#'
#' @template apa_list_arg
#'
#' @param nsim Integer, how many simulation runs will be done.
#'
#' @param save_folder If specified, intermediate simulation results and
#'   processing report will be stored in this folder (defaults to NULL).
#'
#' @param save_simulations If TRUE (default) and if \code{save_folder} is
#'   specified,  individual simulations will be stored (may consume large
#'   storage capacities).
#'
#' @param overwrite If FALSE (default), \code{save_folder} needs to be empty. If
#'   TRUE, all content in \code{save_folder} will be deleted before the new
#'   simulations start.
#'
#' @param buffer_separate Should tree locations in the buffer and cores zones
#'   be randomized separately? Defaults to TRUE.
#'
#' @param parallel If TRUE (default), parallel processing will be used.
#'
#' @param cl If \code{parallel} is TRUE, a Parallel Socket Cluster that was
#'   created with the [parallel::makeCluster()]-function is specified here. If
#'   NULL (default), \code{parallel::makeCluster(no_cores)} will be used.
#'
#' @param no_cores The number of cores that will be used for parallel processing
#'   (default is 1). If NULL is specified, \code{parallel::detectCores() - 1}
#'   will be used to automatically use one core less than available.
#'
#' @param iseed A seed that will be passed to [parallel::clusterSetRNGStream] to
#'   make simulations reproducible when parallel processing is used (default is
#'   42).
#'
#' @details See Glatthorn (2021) for details.
#'
#' @template ref_glatthorn_2021
#'
#' @example examples/example_apa_seg.R
#'
#' @export
apa_seg <- function(apa_list,
                    nsim = 1000, 
                    save_folder = NULL, 
                    save_simulations = TRUE,
                    overwrite = FALSE, 
                    buffer_separate = TRUE, 
                    parallel = TRUE,
                    cl = NULL,
                    no_cores = 1,
                    iseed = 42){
  apa_config <- attr(apa_list, "apa_config")
  apa_properties <- apa_config$apa_properties
  if(!is.null(save_folder)){
    if(!dir.exists(save_folder)){
      dir.create(save_folder)
    }
    previous_runs <- dir(save_folder, full.names = TRUE)
    if(length(dir(save_folder)) != 0 & !overwrite){
      stop("Directory `", save_folder, "` is not empty. Select `overwrite = TRUE`` to delete previous simulation runs.")
    }
    unlink(previous_runs, recursive = TRUE)
  }
  
  sim_idx_nchar <- nchar(nsim)
  sim_idx <- sprintf(paste0("%0", sim_idx_nchar, "d") , 1:nsim)
  if(is.na(apa_properties[1]) || !any(grepl("_ndiv", apa_properties)) ||
     !any(grepl("_pdiv", apa_properties))){
    stop("`ndiv` and `pdiv` are no an apa-properties in `apa_list`. Please use ",
         "`apa_ndiv(apa_list, pdiv = TRUE)` to calculate ndiv before attempting to ",
         "calculate the segregation of ndiv.")
  }
  apa_config <- attr(apa_list, "apa_config")
  apa_properties <- apa_config$apa_properties
  
  drop_tree_dat <- function(x){
    if(!is.null(x$tree_dat)){
      x$tree_dat <- data.frame()
    }
    x
  }
  
  paste0(names(apa_config$dis_trait_column), "_ndiv")
  keep_columns <- c(apa_config$plot_id_column,
                    apa_config$tree_id_column,
                    apa_config$agg_class_column,
                    apa_config$subplot_id_column,
                    paste0(names(apa_config$dis_trait_column), "_ndiv"),
                    paste0(names(apa_config$dis_trait_column), "_pdiv"))
  observed_output <- 
    lapply_deep(apa_list, .f = drop_columns_except, keep_columns = keep_columns,
                .other = drop_tree_dat)
  
  start_time <- Sys.time()
  
  if(parallel){
    if(is.null(cl)){
      if(is.null(no_cores)){
        no_cores <- parallel::detectCores() - 1
      }
      cl <- parallel::makeCluster(no_cores)
    }
    # `loadNamespace("APAtree")` causes some Namespace problems that I completely fail to understand.
    # `library("APAtree)` works fine.
    if(!is.null(iseed)){
      parallel::clusterSetRNGStream(cl, iseed = iseed)
    }
    junk <- parallel::clusterEvalQ(cl, library("APAtree"))
    sim_output <- parallel::parLapply(cl, sim_idx, nsim = nsim, fun = seg_randomize,
                                      apa_list = apa_list, save_folder = save_folder,  save_simulations = save_simulations, 
                                      buffer_separate = buffer_separate, keep_columns = keep_columns,
                                      start_time = start_time)
    parallel::stopCluster(cl)
  }else{
    sim_output <- lapply(sim_idx, nsim = nsim, FUN = seg_randomize,buffer_separate = buffer_separate,
                         apa_list = apa_list, save_folder = save_folder, save_simulations = save_simulations,
                         keep_columns = keep_columns,
                         start_time = start_time)
  }
  if(!is.null(save_folder)){
    temp_files <- grep("temp.txt", dir(save_folder, full.names = TRUE), value = TRUE)
    unlink(temp_files)
  }
  apa_prop_columns <- 
    c(paste0(names(apa_config$dis_trait_column), "_ndiv"),
      paste0(names(apa_config$dis_trait_column), "_pdiv"))
  sim_output <- 
    mapply_deep(c(list(observed_output), sim_output),
                combine_df_cols, 
                value_columns = apa_prop_columns)
  append_to_column_name <- function(x, column_names, append){
    columns_in_x <- intersect(names(x), column_names)
    new_column_names <- paste(columns_in_x, append, sep = "_")
    names(x)[names(x) %in% column_names] <- new_column_names
    x
  }
  sim_output <- lapply_deep(sim_output, .f = seg_orthogonalize)
  sim_output <- 
    lapply_deep(sim_output, .f = append_to_column_name, 
                column_names = apa_prop_columns,
                append = "replications")
  
  if(is.character(save_folder)){
    save(sim_output, file = paste0(save_folder, "/sim_output.rda"))
  }
  seg_columns <- paste0("seg_", names(apa_config$dis_trait_column), "_ndiv")
  apa_list_seg <- 
    mapply_deep(list(apa_list, sim_output), .f = add_col,
                col = seg_columns)
  class(apa_list_seg) <- c("apa_list", "list")
  apa_config_new <- apa_config
  apa_config_new$apa_properties <- c(apa_config_new$apa_properties, seg_columns)
  attr(apa_list_seg, "apa_config") <- apa_config_new
  validate_apa_list(apa_list_seg)
  if(!is.null(save_folder)){
    save(apa_list_seg, file = paste0(save_folder, "/apa_list_seg.rda"))
  }
  apa_list_seg
}

#' @keywords internal
randomize_apa_list <- 
  function(apa_list, 
           buffer_separate = TRUE, 
           save_folder = NULL, 
           suffix = NULL,
           update = TRUE){
    apa_config <- attr(apa_list, "apa_config")
    # Randomization is done over each each element in `plot_id_values`
    if(length(apa_config$plot_id_values) > 1){
      apa_sub_list <- lapply(X = apa_config$plot_id_values, FUN = subset.apa_list,
                             x = apa_list)
      apa_sub_list_sim <- lapply(apa_sub_list, randomize_apa_list, buffer_separate = buffer_separate,
                                 save_folder = NULL, suffix = NULL, update = FALSE)
      apa_list_sim <- mapply_deep(apa_sub_list_sim, .f = function(x){do.call(what = rbind, x)})
      apa_config$randomized <- TRUE
      attr(apa_list_sim, "apa_config") <- apa_config
    }else{
      apa_list_sim <- apa_list
      core_poly <- apa_list_sim$plot_dat[[apa_config$core_column]]
      buffer_poly <- apa_list_sim$plot_dat[[apa_config$buffer_column]]
      if(!buffer_separate){
        core_poly <- buffer_poly
      }
      donut_poly <- core_poly[NULL]
      for(poly_i in seq_along(core_poly)){
        donut_poly_i <- sf::st_difference(buffer_poly[poly_i],
                                          core_poly[poly_i])
        if(length(donut_poly_i) != 0){
          donut_poly <- c(donut_poly, donut_poly_i)
        }
      }
      donut_poly <- sf::st_sfc(donut_poly)
      resample_poly <- c(core_poly, donut_poly)
      if(any(sf::st_overlaps(resample_poly, sparse = FALSE))){
        stop("Some polygons for resampling are overlapping.")
      }
      
      intersect_mat <- sf::st_intersects(apa_list_sim$tree_dat, resample_poly, sparse = FALSE)
      intersect_vec <- col(intersect_mat)
      intersect_vec[!intersect_mat] <- 0 
      intersect_vec <- apply(intersect_vec, 1, max)
      
      resample_freq <- table(intersect_vec)
      if(!"0" %in% names(resample_freq)){
        resample_freq <- c("0" = 0, resample_freq)
      }
      resample_points <- 
        c(sf::st_geometry(apa_list_sim$tree_dat)[intersect_vec == 0],
          sf::st_sample(resample_poly, size = resample_freq[-1]))
      resample_points_sorted <- resample_points[NULL]
      resample_points_sorted[sort(intersect_vec, index.return = TRUE)$ix] <- 
        resample_points
      sf::st_geometry(apa_list_sim$tree_dat) <-  resample_points_sorted
      apa_config$randomized <- TRUE
      attr(apa_list_sim, "apa_config")$randomized <- TRUE
    }
    if(update){
      apa_list_sim <- update_apa(apa_list_sim, force = "all", no_polygon_update = TRUE)
    }
    if(!is.null(save_folder)){
      object_name <- paste0("apa_list_sim_", suffix)
      assign(object_name, apa_list_sim)
      save(list = object_name,
           file = paste0(save_folder, "/sim_", suffix, ".rda"))
    }
    validate_apa_list(apa_list_sim)
  }

#' @keywords internal 
seg_orthogonalize <- 
  function(ndiv, pdiv = NULL) {
    UseMethod("seg_orthogonalize")
  }

#' @keywords internal 
seg_orthogonalize.default <-
  function(ndiv, pdiv){
    ndiv_obs <- ndiv[1]
    pdiv_obs <- pdiv[1]
    if(sum(!is.na(ndiv[-1]) & !is.na(pdiv[-1])) < 3){
      return(NA)
    }
    regression_data <- 
      data.frame(ndiv_sim = ndiv[-1],
                 pdiv_sim = pdiv[-1])
    lm_1 <- stats::lm(ndiv ~ pdiv, data = regression_data)
    pred <- stats::predict(lm_1, newdata = data.frame(pdiv = pdiv_obs), se.fit = TRUE, interval = "prediction")
    ndiv_mean_exp <- pred$fit[1, "fit"]
    ndiv_sig_exp = sqrt(pred$se.fit^2 +pred$residual.scale^2)
    (ndiv_obs - ndiv_mean_exp) / ndiv_sig_exp
  }

#' @keywords internal 
seg_orthogonalize.data.frame <- function(ndiv){
  ndiv_names <- sort(grep("_ndiv", names(ndiv), value = TRUE))
  pdiv_names <- sort(grep("_pdiv", names(ndiv), value = TRUE))
  for(ndiv_i in ndiv_names){
    ndiv_i
    trait_i <- sub("_ndiv", "", ndiv_i)
    pdiv_i <- paste0(trait_i, "_pdiv")
    seg_ndiv_i <- Vectorize(seg_orthogonalize.default)(ndiv[[ndiv_i]], ndiv[[pdiv_i]])
    ndiv[[paste0("seg_", trait_i, "_ndiv")]] <- seg_ndiv_i
  }
  # ndiv <- ndiv[intersect(names(ndiv), c("species_ndiv", "species_pdiv"))]
  # ndiv$seg_ndiv <- seg_ndiv
  ndiv
}

#' @keywords internal 
seg_randomize <- 
  function(sim_i, nsim, apa_list, save_folder = NULL, save_simulations = TRUE,
           buffer_separate = NULL, 
           keep_columns = keep_columns, start_time = NULL){
    {
      if(!is.null(save_folder)){
        proc_file <- file(paste0(save_folder, "/", Sys.getpid() , "_log.txt"),
                          open = "at")
        sink(proc_file, type = "message")
      }
      message("\n----------- SIMULATION ", sim_i, " -----------")
      sim_output_i <- randomize_apa_list(apa_list, buffer_separate = buffer_separate,
                                    save_folder = if(save_simulations){save_folder}else{NULL},
                                    suffix = sim_i)
      if(!is.null(save_folder)){
        sink(type = "message")
        close(proc_file)
        cat(sim_i, "\n", file = paste0(save_folder, "/", Sys.getpid(), "_temp.txt"), append = TRUE)
        temp_files <- grep("temp.txt", dir(save_folder, full.names = TRUE), value = TRUE)
        sim_proc <- as.integer(c(lapply(temp_files, readLines), recursive = TRUE))
        proc_n <- length(sim_proc)
        time_passed <- difftime(Sys.time(), start_time)
        time_total_est <- time_passed / proc_n * nsim
        cat(
          "Starting time of the simulations: ", format(start_time), "\n",
          proc_n, " of ", nsim , " simulations processed (",
          round(100 * proc_n / nsim, 2)," %).\n",
          "Estimated time left: ", round(time_total_est - time_passed, 2), " ",
          attr(time_total_est, "units"), ".\n",
          "Estimated finishing time: ", 
          format(start_time + time_total_est), ".\n",
          sep = "",
          file = paste0(save_folder, "/proc.txt"))
      }
      
      drop_tree_dat <- function(x){
        if(!is.null(x[["tree_dat"]])){
          x$tree_dat <- data.frame()
        }
        x
      }
      sim_output_i <-
        lapply_deep(sim_output_i, .f = drop_columns_except,
                    keep_columns = keep_columns, .other = drop_tree_dat)
    }
    sim_output_i
  }
