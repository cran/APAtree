## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
oldpar <- par(no.readonly = TRUE)

## ----tree data----------------------------------------------------------------
library(APAtree)

# invent some tree data
example_tree_dat <- 
  data.frame(id_tree = as.character(1:9),
             id_plot = "A",
             weight = c(4, 2.6, 4, 5, 3, 5, 3, 4, 5),
             species = c("F", "F", "P", "P", "F", "P", "F", "F", "P"),
             x = c(5, 2, 4, 8, 5, -1, -2, 12, 7),
             y = c(5, 2, 8, 4, 11, 7, 1.5, 7, -1))


## ----tree data formatting-----------------------------------------------------
example_tree_dat <- st_as_sf(example_tree_dat, coords = c("x", "y"))

## ----plot data formatting-----------------------------------------------------
# A test plot area surrounding the threes
plot_poly <- 
  st_polygon(x = list(matrix(c(0, 0, 0, 10, 10, 10, 10, 0, 0, 0),
                             ncol = 2, byrow = TRUE))) %>% 
  st_sfc()

buffer_poly <- 
  st_polygon(x = list(matrix(c(-3, -3, -3, 13, 13, 13, 13, -3, -3, -3),
                             ncol = 2, byrow = TRUE))) %>% 
  st_sfc()

example_plot_dat <- 
  st_sf(id_plot = "A", 
        buffer_geometry = buffer_poly,
        border_geometry = plot_poly)

## ----plot raw data, fig.width = 2, fig.height = 2, fig.align="center"---------
par(mai = rep(0.01, 4))
plot(example_plot_dat$buffer_geometry)
plot(example_plot_dat$border_geometry, add = TRUE)
plot(example_tree_dat$geometry, add = TRUE)

## ----create apa_list----------------------------------------------------------
example_apa_list <- 
  apa_list(tree_dat = example_tree_dat,
           plot_dat = example_plot_dat, 
           plot_id_column = "id_plot",
           tree_id_column = "id_tree",
           core_column = "border_geometry", 
           buffer_column = "buffer_geometry", 
           weight_column = "weight",
           apa_polygon =  TRUE)
print(example_apa_list)

## ----code show APA-maps-------------------------------------------------------
example_apa_list$plot_dat$apa_map

## ----plot apa_list, fig.width = 2, fig.height = 2, fig.align="center"---------
# Default settings
plot(example_apa_list)
# Do not color the critical area
plot(example_apa_list, critical  = NA)
# Use custom colors 
plot(example_apa_list, critical  = NA,
     color_map = 
       data.frame(species = c("P", "F"),
                  color = sf::sf.colors(2, categorical = TRUE)))


## ----apa_size-----------------------------------------------------------------
example_apa_list_2 <- apa_size(example_apa_list)

# The absolute size of the APA-patches within the core plot of all trees (NA means that the tree doesn't have an APA within the core plot).
example_apa_list_2$tree_dat$apa_size

# The relative size of the APA-patches within the core plot of all trees (NA means that the tree doesn't have an APA within the core plot).
example_apa_list_2$tree_dat$apa_size_prop

# Information about the absolute size of the critical area of each tree
example_apa_list_2$tree_dat$critical_area

# Information about the proportion of the total apa size of a tree that is critical
example_apa_list_2$tree_dat$critical

# The critical area may be ignored when calculating the APA-size:
example_apa_list_3 <- apa_size(example_apa_list, edge_correction = "critical")

# In this case, the APA-size gets smaller:
example_apa_list_2$tree_dat$apa_size - example_apa_list_3$tree_dat$apa_size

# The apa_size_total column always gives the overall APA-size irrespective of the used edge-correction method: 
example_apa_list_3$tree_dat$apa_size_total

# The relative apa size is always calculated relative to the APA-size with edge-correction: 
all.equal(example_apa_list_3$tree_dat$apa_size /
            sum(example_apa_list_3$tree_dat$apa_size, na.rm = TRUE),
          example_apa_list_3$tree_dat$apa_size_prop)



## ----apa_ndiv-----------------------------------------------------------------
example_apa_list_4 <- apa_ndiv(example_apa_list, dis_trait_column = "species")

# SpeciesNDiv of the trees (proportion of mixed borders surrounding each tree):
example_apa_list_4$tree_dat$species_ndiv

# SpeciesNDiv at stand-level:
example_apa_list_4$plot_dat$species_ndiv
# Which is the weighted average of all individual SpeciesNDiv values:
weighted.mean(example_apa_list_4$tree_dat$species_ndiv,
              example_apa_list_4$tree_dat$apa_size,
              na.rm = TRUE)
# As the APA-size is required for upscaling, it is always added before calculating NDiv.

## ----simpson------------------------------------------------------------------
# pdiv
(pdiv = example_apa_list_4$plot_dat$species_pdiv)

# Simpson diversity using APA-sizes
species_apa_size <- 
  tapply(example_apa_list_4$tree_dat$apa_size, 
         example_apa_list_4$tree_dat$species,
         sum, na.rm = TRUE)
(simpson <- 1 - sum((species_apa_size / 100)^2))
all.equal(pdiv, simpson)

## ----apa_add_agg_class--------------------------------------------------------
# The aggregation classes are specifies as columns of the tree data
example_apa_list_5 <- 
  apa_add_agg_class(example_apa_list, agg_class_column = "species")

## ----aggregated ndiv----------------------------------------------------------
example_apa_list_6 <-
  apa_ndiv(example_apa_list_5, dis_trait_column = "species")
example_apa_list_6$species[, c("species", "species_ndiv")]

## ---- include = FALSE---------------------------------------------------------
par(oldpar)

