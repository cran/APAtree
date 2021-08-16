# Processing script for Glatthorn et al. 2020 . This script produces all apa
# calculations to calculate the segregation of all plots and
# subplots. WARNING: to process this script may take a while.

library(APAtree)
data(tree_enrico, package = "APAtree")
data(plot_enrico, package = "APAtree")
sf::st_geometry(plot_enrico) <- "border_geometry"


# Create a folder to store results
results_folder <- "./publication_results"
if(!results_folder %in% list.dirs(recursive = TRUE)){
  dir.create(results_folder)
}
sim_folder <- paste0(results_folder, "/simulation_results")
if(!sim_folder %in% list.dirs(recursive = TRUE)){
  dir.create(sim_folder)
}


# run simulations ---------------------------------------------------------

apa_list_enrico <-
  APAtree::apa_list(plot_dat = plot_enrico,
                    tree_dat = tree_enrico,
                    plot_id_column = "id_plot",
                    tree_id_column = "id_tree",
                    core_column = "border_geometry",
                    buffer_column = "buffer_geometry",
                    weight_column = "crown_radius_95",
                    res = .1,
                    apa_polygon = TRUE,
                    apa_properties = c("pdiv"),
                    dis_trait_column = list("species", "height"))

apa_list_enrico_seg <- 
  APAtree::apa_seg(apa_list_enrico, nsim = 1000,
          save_folder = sim_folder,
          overwrite = TRUE,
          parallel = TRUE,
          no_cores = 2, # only using 2 cores, increase to reduce processing time
          iseed = 42)
