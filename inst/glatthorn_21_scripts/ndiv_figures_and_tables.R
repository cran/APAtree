# Processing script for the figures and tables in Glatthorn 2021. All APA
# calculations of simulations are in the script `calculate_simulation_output.R`,
# which has to be processed first. This script will remain static, it will NOT
# be updated to potential future versions of APAtree or to future versions of
# any of the packages APAtree depends on. The sessionInfo() output when the
# script was processed the last time is stored in the file session_info.txt.

library(APAtree)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

results_folder <- "./publication_results"
sim_folder <- paste0(results_folder, "/simulation_results")

# Create a folder to store figures
figure_folder <- paste0(results_folder, "/figures_and_tables")
if(!figure_folder %in% list.dirs(recursive = TRUE)){
  dir.create(figure_folder)
}

# store old par settings

oldpar <- par(no.readonly = TRUE)

# data loading ---------------------------------------------

# package data
data(tree_enrico, package = "APAtree")
data(plot_enrico, package = "APAtree")
data(subplot_enrico, package = "APAtree")

# simulation data calculated in `calculate_simulation_output.R`
load(paste0(sim_folder, "/sim_output.rda"))
load(paste0(sim_folder, "/apa_list_seg.rda"))

# Data formatting --------------

sim_stand <- 
  sim_output$plot_dat %>% 
  rowwise() %>% 
  do({tibble(id_plot = .$id_plot,
             height_ndiv = c(.$height_ndiv_replications, recursive = TRUE),
             species_ndiv = c(.$species_ndiv_replications, recursive = TRUE),
             RD = c(.$height_pdiv_replications, recursive = TRUE),
             SD = c(.$species_pdiv_replications, recursive = TRUE),
  )}) %>% 
  ungroup() %>% 
  group_by(id_plot) %>% 
  mutate(n_sim = 1:n())


color_map <- 
  data.frame(species = c("Fagus sylvatica", "Pseudotsuga menziesii"),
             species_color = c("#7d5831", "#bcc746"))

color_map <- 
  unique(select(st_drop_geometry(apa_list_seg$tree_dat), species)) %>% 
  mutate(color = case_when(species == "Fagus sylvatica" ~ "#7d5831",
                           species == "Pseudotsuga menziesii" ~ "#bcc746",
                           species == "Picea abies" ~ "#b6d7a8",
                           TRUE ~ "#c5d9ed"))
plot_limits <-
  plot_enrico %>% 
  split(.$id_plot) %>% 
  map(st_bbox)

# Figure 1 ----------------------------------------------------------------

plot_limits[["8.2"]] <-
  plot_enrico %>%
  filter(id_plot == "8.2") %>%
  st_bbox
plot_limits[["5.2"]] <-
  plot_enrico %>%
  filter(id_plot == "5.2") %>%
  st_bbox


{
  svg(file = paste0(figure_folder , "/figure_1.svg"),
      width = 140/25.4,
      height = 2.6,
      pointsize = 12,
      bg = NA)
  par(mfrow = c(1, 2))
  par(omi = c(0, 0, 0, 0.015))
  par(mai = c(0.01, .17, 0.01, 0.01))
  plot.new()
  plot.window(xlim = plot_limits[["5.2"]][c(1, 3)], ylim = plot_limits[["5.2"]][c(2, 4)],xaxs = "i", yaxs = "i", asp = 1)
  plot(apa_list_seg, "5.2",
       color_map = color_map, add_legend = FALSE,
       tree_size_column = "dbh", add_plot_id_values = FALSE, add = TRUE)
  mtext("A", side = 3, adj = -0.07, font = 2, line = -1.1)
  plot.new()
  par(mai = c(0.01, .02, 0.01, 0.01))
  plot.window(xlim = plot_limits[["8.2"]][c(1, 3)], ylim = plot_limits[["8.2"]][c(2, 4)], xaxs = "i", yaxs = "i", asp = 1)
  plot(apa_list_seg, "8.2",
       color_map = color_map, add_legend = FALSE,
       tree_size_column = "dbh", add_plot_id_values = FALSE, add = TRUE)
  mtext("B", side = 3, adj = 0.05, font = 2, line = -1.1)
  dev.off()
}


# Figure S1 ---------------------------------------------------------------

{
  tiff(file = paste0(figure_folder , "/figure_S1.tiff"),
       width = 175/25.4,
       height = 175/25.4,
       pointsize = 12,
       bg = NA,
       units = "in",
       res = 150,
       compression = "lzw")
  par(mfrow = c(4, 4))
  par(omi = c(0, 0, 0, 0))
  for(i in seq_along(plot_limits)){
    print(i)
    par(mai = rep(0, 4))
    plot.new()
    par(mai = c(0.01, .01, 0.13, 0.01))
    plot.window(xlim = plot_limits[[i]][c(1, 3)],
                ylim = plot_limits[[i]][c(2, 4)],
                xaxs = "i", yaxs = "i", asp = 1)
    mtext(names(plot_limits)[i], side = 3, adj = 0.6, font = 2, line = 0, at = 0, cex = .8)
    plot(apa_list_seg, i,
         color_map = color_map, add_legend = FALSE,
         tree_size_column = "dbh", add_plot_id_values = FALSE, add = TRUE)
  }
  dev.off()
}

# Figure 2 ----------------------------------------------------------------

# Examplary dataset with two trees
example_tree_dat <- 
  tibble(id_tree = as.character(1:9),
         id_plot = "A",
         competitiveness = c(4, 2.6, 4, 5, 3, 5, 3, 4, 5),
         position = c(rep("plot", 4), rep("buffer", each = 5)),
         x = c(5, 2, 4, 8, 5, -1, -2, 12, 7),
         y = c(5, 2, 8, 4, 11, 7, 1.5, 7, -1)) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  rename(stem_position = geometry)

# A virtual plot area surrounding both threes
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



# apa_list - 1 m resolution to visualize raster
example_apa_list <-
  apa_list(plot_dat = example_plot_dat, plot_id_column = "id_plot",
           tree_dat = filter(example_tree_dat, position == "plot"),
           core_column = "border_geometry",
           weight_column = "competitiveness", res = 1, tree_id_column = "id_tree",
           apa_polygon = TRUE)

# centers of the grid cells
grid_cells <-
  st_as_sf(rasterToPoints(example_apa_list$plot_dat$apa_map[[1]], spatial = TRUE))$geometry

{
  svg(filename = paste0(figure_folder , "/figure_2_R.svg"),
      width = 90/25.4,
      height = 1.7,
      pointsize = 12,
      bg = NA)
  par(oma = rep(0.01, 4))
  par(mai = c(0.01, 0.12, 0.01, 0.01))
  layout(rbind(1:2))
  xlim = c(-0.1, 10.1)
  ylim = c(-0.1, 10.1)
  
  # Panel A
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1, xaxs = "i", yaxs = "i")
  plot(plot_poly, add = TRUE, lwd = 6, border = gray(.7), xpd = NA)
  plot(grid_cells, add = TRUE, pch = 16, cex = .3)
  plot(example_apa_list$tree_dat$apa_polygon, add = TRUE, lwd = 2)
  example_tree_dat %>% 
    filter(position == "plot") %>% 
    {symbols(x = st_coordinates(.$stem_position), add = TRUE, pch = 16,
             circles = .$competitiveness/10, bg = color_map$color[c(1, 1, 2, 2)], inches = FALSE)}
  mtext("A", side = 3, adj = -0.08, font = 2, line = -0.9)
  
  # Panel B
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1, yaxs = "i", xaxs = "i")
  plot(plot_poly, add = TRUE, lwd = 6, border = gray(.7), xpd = NA)
  plot(example_apa_list$tree_dat$apa_polygon, add = TRUE, lwd = 2,
       col = paste0(substr(color_map$color[c(2, 2, 1, 1)], 0, 7), "44"))
  example_tree_dat %>% 
    filter(position == "plot") %>% 
    {symbols(x = st_coordinates(.$stem_position), add = TRUE, pch = 16,
             circles = .$competitiveness/10, bg = color_map$color[c(2, 2, 1, 1)], inches = FALSE)}
  mtext("B", side = 3, adj = -0.08, font = 2, line = -.9)
  
  boundaries <-
    st_intersection(example_apa_list$tree_dat$apa_polygon[1],
                    example_apa_list$tree_dat$apa_polygon[2:4])
  boundaries <-
    st_sf(focal_tree = 1, neighbor_tree = 2:4, boundaries)
  plot(boundaries$boundaries[1:3], col = color_map$color[2], add = TRUE, lwd = 4)
  plot(boundaries$boundaries[2:3], col = color_map$color[1], lty = 3, add = TRUE, lwd = 4)
  
  dev.off()
}

# The Figure as plotted by R was modified with inkscape for the publication

# Figure 4 ----------------------------------------------------------------

sim_stand_42 <- filter(sim_stand, id_plot == "4.2")
{
  svg(paste0(figure_folder, "/figure_4_R.svg"), 
      width = 90/25.4, height = 5, pointsize = 12)
  par(mai = c(.6, .6, .01, .01))
  par(mgp = c(2, 1, 0))
  plot(species_ndiv ~ SD,
       data = filter(sim_stand_42, n_sim != 1),
       pch = 16, col = gray(.85),
       xlab = "Simpson diversity", 
       yaxt = "n", xaxt = "n",
       ylab = expression("SpeciesNDiv"["stand"]),
       asp = 1, cex = .5,
       ylim = c(0.325, 0.66),
       xlim = range(sim_stand_42$SD),
       bty = "n", xaxs = "i", yaxs = "i")
  axis(1, at = c(0.4, 0.5, 0.6, .7))
  axis(2, at = c(0.30, 0.4, .5, .6, .7))
  
  sim_beech_lm <- 
    lm(species_ndiv ~ SD,
       data = filter(sim_stand_42, n_sim != 1))
  
  species_ndiv_pred <- 
    sim_stand_42 %>% 
    {tibble(SD = seq(min(.$SD), max(.$SD), length.out = 100))}
  
  species_ndiv_pred <- 
    cbind(species_ndiv_pred, 
          predict(sim_beech_lm, newdata = species_ndiv_pred, se.fit = TRUE, interval = "prediction")[c(1, 2, 4)])
  species_ndiv_pred <-
    species_ndiv_pred %>%
    mutate(
      cond_sd = sqrt(se.fit^2 + residual.scale^2)
    )
  
  species_ndiv_pred %>% 
    {polygon(c(.$SD, rev(.$SD)),
             y =c(.$fit.fit, rev(.$fit.fit)) + c(.$cond_sd, rev(-1*.$cond_sd)),
             col = "#228B2244", border = NA)}
  lines(fit.fit ~ SD, data = species_ndiv_pred, lty = 1, lwd = 2)
  
  sim_stand_42 %>% 
    {points(.$species_ndiv[1] ~ .$SD[1],
            pch = 16, cex = 1.5)}
  dev.off()
}

# Figure 4 caption reference --------------------------------------------------------

apa_list_seg$plot_dat %>% 
  filter(id_plot == "4.2") %>% 
  pluck("seg_species_ndiv") %>% 
  round(2)

# Figure 5 ----------------------------------------------------------------

sim_stand <-  
  sim_stand %>% 
  ungroup() %>% 
  mutate(doug_spruce = 3 - as.numeric(substr(.$id_plot, 3, 3))/2)

{
  tiff(filename = paste0(figure_folder , "/figure_5.tiff"),
       width = 140,
       height = 70,
       pointsize = 12,
       bg = NA,
       res = 1000,
       compression = "lzw",
       units = "mm")
  par(mfrow = c(1, 2))
  par(omi = rep(0.01, 4))
  par(mai = c(.5, .5, .02, .01))
  par(mgp = c(1.5, .5, 0))
  par(tcl = -.4)
  sim_stand %>% 
    filter(n_sim != 1) %>% 
    {plot(species_ndiv ~ SD , data = ., pch = 16,
          #col = paste0(color_map$species_color, "05")[doug_spruce],
          col = gray(0.6, .04),
          xlab = "Simpson diversity", ylab = expression(SpeciesNDiv["stand"]))}
  mtext("A", side = 3, adj = 0.05, font = 2, line = -1.3)
  
  lm_1 <- 
    sim_stand %>% 
    filter(n_sim != 1) %>% 
    group_by(id_plot, doug_spruce) %>% 
    summarise_at(c("species_ndiv", "SD"),
                 mean) %>% 
    {lm(species_ndiv ~ SD, 
        data = .)}
  
  lm_1_pred <-
    sim_stand %>%
    ungroup() %>%
    pluck("SD") %>%
    range() %>%
    tibble(SD = .) %>%
    {add_column(., fit = predict(lm_1, .))}
  lines(x = lm_1_pred$SD, y = lm_1_pred$fit, lwd = 3, lty = 2,
        #col = gray(0.8)
        col = "black")
  
  lm_2 <- 
    sim_stand %>%
    ungroup() %>%
    lm(formula = species_ndiv ~ SD * id_plot)
  
  lm2_pred <- 
    sim_stand %>% 
    group_by(id_plot) %>% 
    filter(n_sim != 1) %>% 
    do(
      tibble(
        min_max = c("min", "max"),
        SD = c(min(.$SD), max(.$SD)),
        avg = mean(.$SD)))
  lm2_pred <- 
    lm2_pred %>% 
    {add_column(., fit = predict(lm_2, newdata = .))}
  lm2_pred %>% 
    split(.$id_plot) %>% 
    map(~lines(.x$SD, .x$fit))
  
  sim_stand %>% 
    filter(n_sim == 1) %>% 
    {points(.$species_ndiv ~ .$SD,
            #data = m_prop_apa_beech_core,
            pch = 21,
            bg = gray(0.6))}
  par(mai = c(.5, .55, .02, .01))
  apa_list_seg$plot_dat %>% 
    {plot(seg_species_ndiv ~ species_pdiv , data = .,
          pch = 21, bg = gray(0.6),
          xlab = "Simpson diversity", ylab = expression(seg(SpeciesNDiv["stand"])))}
  mtext("B", side = 3, adj = 0.05, font = 2, line = -1.3)
  
  
  dev.off()
  }

# Figure 5 caption and text references -----------------------------------------------

lm_SD <- 
  sim_stand %>% 
  filter(n_sim != 1) %>% 
  group_by(id_plot, doug_spruce) %>% 
  summarise_at(c("species_ndiv", "SD"),
               mean) %>% 
  {lm(species_ndiv ~ SD, 
      data = .)}

# Simulation model coefficients (caption)
round(coef(lm_SD), 2)

# Difference between simulated NDiv and SD values + CI (results section):
diff_sim_NDiv_SD <- t.test(lm_SD$model$species_ndiv - lm_SD$model$SD)
round(diff_sim_NDiv_SD$estimate, 3)
round(diff_sim_NDiv_SD$conf.int, 3)

# Simulation model correlation coefficient (caption and results)
round(sqrt(summary(lm_SD)$r.squared), 3)

#Correlation between SpeciesNDiv and Simpson diversity (SpeciesPDiv) of the
#actual observations (caption and results):
sim_stand %>% 
  filter(n_sim == 1) %>% 
  {cor(.$species_ndiv, .$SD)} %>% 
  round(2)

# Correlation between seg_SpeciesNDiv and Simpson diversity (SpeciesPDiv)
# (caption and results):
apa_list_seg$plot_dat %>% 
  {cor(.$seg_species_ndiv, .$species_pdiv)} %>% 
  round(2)

# Figure 6 ----------------------------------------------------------------

sim_stand <- 
  sim_stand %>% 
  ungroup() %>% 
  mutate(doug_spruce = as.numeric(substr(.$id_plot, 3, 3))/2)

{
  tiff(filename = paste0(figure_folder , "/figure_6.tiff"),
       width = 140,
       height = 70,
       pointsize = 12,
       bg = NA,
       res = 1000,
       compression = "lzw",
       units = "mm")
  par(mfrow = c(1, 2))
  par(omi = rep(0.01, 4))
  par(mai = c(.5, .5, .02, .01))
  par(mgp = c(1.5, .5, 0))
  par(tcl = -.4)
  sim_stand %>% 
    filter(n_sim != 1) %>% 
    {plot(height_ndiv ~ RD , data = ., pch = 16,
          #col = ccolor_map$species_color[3 - doug_spruce],
          col = gray(0.6, .04),
          xlab = "Rao diversity", ylab = expression(HeightNDiv["stand"]))}
  mtext("A", side = 3, adj = 0.05, font = 2, line = -1.3)
  
  lm_1 <- 
    sim_stand %>% 
    filter(n_sim != 1) %>% 
    group_by(id_plot, doug_spruce) %>% 
    summarise_at(c("height_ndiv", "RD"),
                 mean) %>% 
    {lm(height_ndiv ~ RD, 
        data = .)}
  summary(lm_1)
  
  lm_1_pred <-
    sim_stand %>%
    ungroup() %>%
    pluck("RD") %>%
    range() %>%
    tibble(RD = .) %>%
    {add_column(., fit = predict(lm_1, .))}
  lines(x = lm_1_pred$RD, y = lm_1_pred$fit, lwd = 3, lty = 2,
        #col = gray(0.8)
        col = "black")
  
  lm_2 <- 
    sim_stand %>%
    ungroup() %>%
    lm(formula = height_ndiv ~ RD * id_plot)
  
  lm2_pred <- 
    sim_stand %>% 
    group_by(id_plot) %>% 
    filter(n_sim != 1) %>% 
    do(
      tibble(
        min_max = c("min", "max"),
        RD = c(min(.$RD), max(.$RD)),
        avg = mean(.$RD)))
  lm2_pred <- 
    lm2_pred %>% 
    {add_column(., fit = predict(lm_2, newdata = .))}
  lm2_pred %>% 
    split(.$id_plot) %>% 
    map(~lines(.x$RD, .x$fit))
  
  sim_stand %>% 
    filter(n_sim == 1) %>% 
    {points(.$height_ndiv ~ .$RD,
            #data = m_prop_apa_beech_core,
            pch = 21, bg = gray(0.6))}
  par(mai = c(.5, .55, .02, .01))
  apa_list_seg$plot %>% 
    {plot(seg_height_ndiv ~ height_pdiv , data = .,
          pch = 21, bg = gray(.6),
          xlab = "Rao diversity", ylab = expression(seg(HeightNDiv["stand"])))}
  mtext("B", side = 3, adj = 0.05, font = 2, line = -1.3)
  
  apa_list_seg$plot %>% 
    {lm(seg_height_ndiv ~ height_pdiv, data = .)} %>% 
    summary()
  
  
  dev.off()
  }

# Figure 6 caption and text references -----------------------------------------------

lm_RD <- 
  sim_stand %>% 
  filter(n_sim != 1) %>% 
  group_by(id_plot, doug_spruce) %>% 
  summarise_at(c("height_ndiv", "RD"),
               mean) %>% 
  {lm(height_ndiv ~ RD, 
      data = .)}

# Simulation model coefficients (caption)
round(coef(lm_RD), 3)

# Difference between simulated NDiv and RD values + CI (results section):
diff_sim_NDiv_RD <- t.test(lm_RD$model$height_ndiv - lm_RD$model$RD)
round(diff_sim_NDiv_RD$estimate, 3)
round(diff_sim_NDiv_RD$conf.int, 3)

# Simulation model correlation coefficient (caption)
round(sqrt(summary(lm_RD)$r.squared), 3)

#Correlation between SpeciesNDiv and Rao diversity (SpeciesPDiv) of the
#actual observations (caption):
sim_stand %>% 
  filter(n_sim == 1) %>% 
  {cor(.$height_ndiv, .$RD)} %>% 
  round(2)

#Correlation between seg_HeightNDiv and Rao diversity (SpeciesPDiv)
#(caption and text):
apa_list_seg$plot_dat %>% 
  {cor(.$seg_height_ndiv, .$height_pdiv)} %>% 
  round(2)

# restore old par settings
par(oldpar)

# Table 1 -------------------------------------------------------

plot_area_tbl <- 
  subplot_enrico %>%
  filter(object == "plot") %>%
  transmute(id_plot, plot_area = st_area(geometry)) %>% 
  mutate(plot_area = units::set_units(plot_area, "m^2"),
         plot_area = units::set_units(plot_area, "ha")) %>% 
  st_set_geometry(NULL)

basic_mixture_properties <- 
  tree_enrico %>% 
  st_set_geometry(NULL) %>% 
  left_join(plot_area_tbl, by = "id_plot") %>% 
  mutate(dbh = units::as_units(dbh, "cm")) %>% 
  filter(position == "plot") %>% 
  group_by(id_plot) %>% 
  summarize(
    avg_dbh = mean(dbh),
    sd_dbh = sd(dbh),
    avg_height = mean(height),
    sd_height = sd(height),
    tree_n = n() / plot_area[1],
    ba = units::set_units(sum(.25 * pi * dbh ^ 2), "m^2") / plot_area[1],
    beech_prop_ba = sum(.25 * pi * dbh[species == "Fagus sylvatica"] ^ 2) /plot_area[1] / ba)

apa_prop_stand <- 
  apa_list_seg$plot_dat %>% 
  st_drop_geometry() %>%
  select(id_plot, species_ndiv, seg_species_ndiv, height_ndiv, seg_height_ndiv)

apa_list_enrico <-
  APAtree::apa_list(plot_dat = plot_enrico,
                    tree_dat = tree_enrico,
                    plot_id_column = "id_plot",
                    tree_id_column = "id_tree",
                    weight_column = "crown_radius_95",
                    res = .1,
                    agg_class_column = "species",
                    subplot_dat = list(subplot = subplot_enrico),
                    subplot_id_column = c(subplot = "id_subplot"),
                    apa_polygon = FALSE,
                    apa_properties = "apa_size")
apa_prop_beech <- 
  apa_list_enrico$subplot_dat$subplot$species %>% 
  filter(grepl("plot", id_subplot), species == "Fagus sylvatica") %>% 
  select(id_plot, beech_prop_apa = apa_size_prop)

table_1_dat <- 
  basic_mixture_properties %>% 
  left_join(apa_prop_beech) %>% 
  left_join(apa_prop_stand) 

table_1_formatted <- 
  table_1_dat %>% 
  mutate_at(2:5, sprintf, fmt = "%.1f") %>% 
  mutate(tree_n = sprintf(fmt = "%.0f", tree_n)) %>% 
  mutate_if(is.numeric, sprintf, fmt = "%.2f") %>% 
  unite(col = "dbh", avg_dbh, sd_dbh, sep = " (± ", remove = TRUE) %>% 
  unite(col = "height", avg_height, sd_height, sep = " (± ", remove = TRUE) %>% 
  mutate_at(c("dbh", "height"), paste0, ")")

writexl::write_xlsx(table_1_formatted, paste0(figure_folder, "/table_1_formatted.xlsx"))

