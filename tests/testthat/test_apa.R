library(APAtree)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(testthat)

# Create test dataset with some trees ------------
test_trees <-
  data.frame(id_tree = as.character(1:9),
             id_plot = "A",
             crown_radius_95 = c(4.02, 2.6, 4.01, 5.01, 3.01, 5, 3.02, 4.03, 5),
             # note: if a grid cell center falls right in between two trees with
             # equal weight, the selection of the tree is platform dependent, I
             # assume due to small differences of floating point operations. To
             # avoid this during testing, each tree gets a slightly different
             # weight.
             position = c(rep("plot", 4), rep("buffer", each = 5)),
             species = c("F", "F", "P", "P", "F", "P", "F", "F", "P"),
             x = c(5, 2, 4, 8, 5, -1, -2, 12, 7),
             y = c(5, 2, 8, 4, 11, 7, 1.5, 7, -1)) %>%
  st_as_sf(coords = c("x", "y")) %>%
  rename(stem_position = geometry)

# A test plot area surrounding the threes
plot_poly <-
  st_polygon(x = list(matrix(c(0, 0, 0, 10, 10, 10, 10, 0, 0, 0),
                             ncol = 2, byrow = TRUE))) %>%
  st_sfc()

buffer_poly <-
  st_polygon(x = list(matrix(c(-3, -3, -3, 13, 13, 13, 13, -3, -3, -3),
                             ncol = 2, byrow = TRUE))) %>%
  st_sfc()

test_plot <-
  st_sf(id_plot = "A",
        buffer_geometry = buffer_poly,
        border_geometry = plot_poly)

# subplot data for testing

# polygon subplots
test_subplot <-
  st_sf(id_subplot = "A1",
        id_plot = "A",
        subplot_geometry = st_buffer(test_plot$border_geometry, -2))
# locations for neighborhood analyses
test_locations <-
  st_sf(id_location = "A2",
        id_plot = "A",
        location_geometry = st_sfc(st_point(x = c(5, 5))))
test_subplot_dat <- list(test_subplot = test_subplot,
                         test_locations = test_locations)

# build an apa_list
test_apa_list <-
  apa_list(tree_dat = test_trees,
           plot_dat = test_plot,
           plot_id_column = "id_plot",
           tree_id_column = "id_tree",
           core_column = "border_geometry",
           buffer_column = "buffer_geometry",
           weight_column = "crown_radius_95",
           apa_polygon = TRUE)

test_apa_list <-
  apa_add_agg_class(test_apa_list,
                    agg_class_column = "species",
                    apa_polygon = TRUE)

test_apa_list <-
  apa_add_subplot_dat(test_apa_list,
                      subplot_dat = test_subplot_dat,
                      subplot_id_column = c(test_subplot = "id_subplot",
                                            test_locations = "id_location"),
                      radius = 4)

# test apa_size output ------------------

test_that("apa_size output is correct", {
  test_apa_size <-
    apa_size(test_apa_list, edge_correction = "critical")
  expect_equal(test_apa_size$tree_dat$apa_size_total,
               c(16, 13, 17, 25, 2, 13, NA, 4, 10))
  expect_equal(test_apa_size$tree_dat$apa_size,
               c(16, 12, 16, 24, 1, 11, NA, 1, 9))
  expect_equal(test_apa_size$tree_dat$apa_size_prop,
               c(0.177777777777778, 0.133333333333333, 0.177777777777778, 0.266666666666667,
                 0.0111111111111111, 0.122222222222222, NA, 0.0111111111111111,
                 0.1))
  expect_equal(test_apa_size$plot_dat$apa_size_total, 100)
  expect_equal(test_apa_size$plot_dat$apa_size, 90)
  expect_equal(test_apa_size$species$apa_size, c(30, 60))
  expect_equal(test_apa_size$species$apa_size_total, c(35, 65))
  expect_equal(test_apa_size$subplot_dat$test_subplot$tree_dat$apa_size,
               c(16, 3, 7, 10))
  expect_equal(test_apa_size$subplot_dat$test_locations$tree_dat$apa_size,
               c(16, 5, 11, 14, 3, 3))
})

# test apa_size output ------------------

test_that("apa_size output is correct", {
  test_apa_ndiv <-
    apa_ndiv(test_apa_list,
             dis_trait_column = "species",
             pdiv = TRUE,
             edge_correction = "critical")
  expect_equal(test_apa_ndiv$tree_dat$species_ndiv,
               c(0.8, 0.333333333333333, 0.615384615384615, 0.529411764705882,
                 1, 0.5, NA, 1, 0.25))
  expect_equal(test_apa_ndiv$tree_dat$species_pdiv,
               c(2/3, 2/3, 1/3, 1/3, 2/3, 1/3, NA, 2/3, 1/3))
  expect_equal(test_apa_ndiv$plot_dat$species_ndiv, 0.545578179989945)
  expect_equal(test_apa_ndiv$plot_dat$species_pdiv, 0.444444444444444)
  expect_equal(test_apa_ndiv$species$species_ndiv, c(0.626666666666667, 0.505033936651584))
  expect_equal(test_apa_ndiv$species$species_pdiv, c(2/3, 1/3))
  expect_equal(test_apa_ndiv$subplot_dat$test_subplot$tree_dat$species_ndiv,
               c(0.764705882352941, 0, 0.857142857142857, 0.875))
  expect_equal(test_apa_ndiv$subplot_dat$test_subplot$test_subplot$species_ndiv,
               0.749591503267974)
  expect_equal(test_apa_ndiv$subplot_dat$test_locations$tree_dat$species_ndiv,
               c(0.8, 0.6, 0.611111111111111, 0.55, 0.555555555555556, 0.4))
  expect_equal(test_apa_ndiv$subplot_dat$test_locations$test_locations$species_ndiv,
               0.6363248)
})
