#'Tree-level data of 16 forest stands in Lower Saxony, Germany
#'
#'This dataset contains tree-level data from the Enrico-project from 16 mixed
#'forest stands in Lower Saxony, Germany.
#'
#'@format A \code{data.frame} of 2813 trees on 16 plots. Tree coordinates are
#'  stored as 'simple feature column' (see the [`sf`][sf::sf]-package for more
#'  details):
#'
#'  \describe{
#'
#'  \item{id_tree}{Unique tree id.}
#'
#'  \item{id_plot}{Plot id's}
#'
#'  \item{species}{Name of the tree species}
#'
#'  \item{position}{Does the tree stand in the core plot ('plot') or in the 10 m
#'  buffer zone around the core plot ('buffer').}
#'
#'  \item{dbh}{Diameter at breast height (1.3 m) in cm.}
#'
#'  \item{height}{Top height of the tree in m. Estimated with a stand height
#'  curve.}
#'
#'  \item{crown_radius_95}{Expected 95 % quantile of the crown radius of a tree
#'  of this species and diameter. Estimated according to allometric equations
#'  (Pretzsch et al. 2015)}
#'
#'  \item{crown_radius_95}{Expected median of the crown radius of a tree of this
#'  species and diameter. Estimated according to allometric equations (Pretzsch
#'  et al. 2015)}
#'
#'  \item{tree_geometry}{A \code{sfc_POINT}-column specifying tree coordinates
#'  relative to the plot center in meter.}
#'
#'  }
#' @template ref_ammer_2020
#'
#' @template ref_pretzsch_2015
#'
#' 
#'@md
"tree_enrico"



#' Plot-level data of 16 forest stands in Lower Saxony, Germany
#'
#' This dataset contains plot-level data from the Enrico-project from 16 mixed
#' forest stands in Lower Saxony, Germany.
#'
#' @format A \code{data.frame} of 16 plots. Plot coordinates are stored as
#'   'simple feature column' (see the [`sf`][sf::sf]-package for more details):
#'
#'   \describe{
#'
#'   \item{id_plot}{Unique plot id's.}
#'
#'   \item{border_geometry}{A \code{sfc_POLYGON}-column specifying the
#'   coordinates of the plot border relative to the plot center.}
#'
#'   \item{buffer_geometry}{A \code{sfc_POLYGON}-column specifying the
#'   coordinates of the plot buffer (10 m distance from \code{boder_geometry})
#'   relative to the plot center.}
#'   }
#'
#' @template ref_ammer_2020
#'
#'@md
"plot_enrico"


#' Subplot-level data of 16 forest stands in Lower Saxony, Germany
#'
#' This dataset contains subplot-level data from the Enrico-project from 16
#' mixed forest stands in Lower Saxony, Germany.
#'
#' @format A \code{data.frame} of 352 subplots on 16 plots. Plot coordinates are
#'   stored as 'simple feature column' (see the [`sf`][sf::sf]-package for more
#'   details):
#'
#'   \describe{
#'
#'   \item{id_subplot}{Unique subplot id's.}
#'
#'   \item{id_plot}{Plot id's.}
#'
#'   \item{object}{Name of the subplots ('plot', 'buffer', 'quadrant',
#'   'sixteenth').}
#'
#'   \item{geometry}{A \code{sfc_POLYGON}-column specifying the coordinates of
#'   the subplots relative to the plot center.}
#'
#'   }
#'
#' @template ref_ammer_2020
#'
#'@md
"subplot_enrico"
