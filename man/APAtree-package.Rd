\name{APAtree-package}
\alias{APAtree-package}
\alias{APAtree}
\docType{package}
\title{
  Generation and analysis of maps of the area potentially available to trees (APA-maps) of forest stands
}
\description{
  APA-maps are calculated from tree inventory data of forest stands. Tree coordinates, the competitive ability of trees and information about plot borders are used as inputs to derive APA-maps using multiplicatively weighted Voronoi diagrams. Functions to calculate important characteristics from APA-maps are provided (e.g., size of APA-patches of individual trees and tree species, neighborhood diversity).
}
\details{
  The core function of the \code{APAtree}-package is the \code{\link{apa_list}} function that creates APA-maps from tree inventory data and stores it together with the input data and APA-properties in an object of class \code{apa_list}.
  
  For a more detailed description of the package see the vignette (\code{vignette("APAtree-vignette", package = "APAtree")})
  
  Additionally, the scripts to produce some of the Figures in Glatthorn (accepted) can be found in the folder \code{paste0(path.package("APAtree"), "/glatthorn_2021")}.
}
\author{
Jonas Glatthorn, jonas.glatthorn@posteo.de.

Maintainer: Jonas Glatthorn <jonas.glatthorn@posteo.de>
}
\references{
{ Glatthorn, Jonas (accepted). A spatially explicit index for
tree species or trait diversity at neighborhood and stand level. Ecological Indicators.}

{Gspaltl, M., Sterba, H., & O'hara, K. L. (2012). The relationship
between available area efficiency and area exploitation index in an even-aged
coast redwood (Sequoia sempervirens) stand. Forestry, 85(5), 567-577.}

{Römisch,K. (1995) Durchmesserwachstum und ebene Bestandesstruktur am Beispiel der Kiefernversuchsfläche Markersbach. In Deutscher Verband forstl.
Forschungsanstalten, Sektion Biometrie und Informatik. Gottfried Hempel (ed.) Tagung Tharanth/Grillenburg, Vol. 8, pp. 84–103.}


}
\keyword{ package }
