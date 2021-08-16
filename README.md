# APAtree

The `area potentially available` of trees is calculated from
  mapped forest stands. This is done by computing a rasterized version of
  'Weighted Voronoi Diagrams' using a an approximation of the trees competitive
  ability (e.g., crown radius, leaf area) as weight. The main output are Raster*
  objects from the raster package that are stored together with the raw data in
  apa_list's, the main class of the APAtree package. Aggregation functions are
  provided to calculate the most important stand characteristics based on APA-maps.
  
## Install

Currently, only installation from github is possible
```r
library(devtools)
install_github("s-spatial/sf")
```

## Documentation

After Installing, see the package Vignette for details (`vignette("APAtree-vignette", package = "APAtree")`)

Additionally, the scripts to produce some of the Figures in Glatthorn (accepted) can be found
in the folder `paste0(path.package("APAtree"), "/glatthorn_2021")`.

## References

Glatthorn, Jonas (accepted): A spatially explicit index for tree species or trait diversity
at neighborhood and stand level. Ecological Indicators.

Gspaltl, M., Sterba, H., & O’hara, K. L. (2012). The relationship
between available area efficiency and area exploitation index in an even-aged
coast redwood (Sequoia sempervirens) stand. Forestry, 85(5), 567-577.

Römisch,K. (1995) Durchmesserwachstum und ebene Bestandesstruktur am Beispiel der Kiefernversuchsfläche Markersbach. In Deutscher Verband forstl.
Forschungsanstalten, Sektion Biometrie und Informatik. Gottfried Hempel (ed.) Tagung Tharanth/Grillenburg, Vol. 8, pp. 84–103.

