#' Add fixed Kernel Utilisation Distribution estimates to 3D rayshader plot
#'
#' @description Helper function to add 3DKUD a 3D rayshader plot
#'
#' @param ras raster object with bathymetry. This is the raster used to create the hillshade nad heightmap of the rayshader 3D plot
#'   tag detection data, metadata and station information
#' @param det ...
#' @param zscale scaling factor for the depth, indicating depth exaggeration
#' @param lonlat ...
#' @param mul ...
#' @param ... additional arguments for \code{\link{plot.kde}} function
#'
#' @return Adds 3DKUD voxels to an existing 3D rayshader plot
#'
#' @import ks
#' @import rayshader
#' @import raster
#' @importFrom raster extent
#' @importFrom raster pointDistance
#' @importFrom rayshader %>%
#' @importFrom ks Hpi
#' @importFrom ks kde
#' @importFrom ks plot.kde
#'
#' @export add_fkud
#'
#' @examples
#' TBD
#'


add_fkud <- function(ras,
                     det,
                     zscale,
                     lonlat = FALSE,
                     mul = 1,
                     ...) {
  e <- raster::extent(ras)
  cell_size_x <-
    raster::pointDistance(c(e@xmin, e@ymin), c(e@xmax, e@ymin), lonlat = lonlat) / ncol(ras)
  cell_size_y <-
    raster::pointDistance(c(e@xmin, e@ymin), c(e@xmin, e@ymax), lonlat = lonlat) / nrow(ras)
  distances_x <-
    raster::pointDistance(c(e@xmin, e@ymin), cbind(det$lon, rep(e@ymin, nrow(det))), lonlat = lonlat) / cell_size_x
  distances_y <-
    raster::pointDistance(c(e@xmin, e@ymin), cbind(rep(e@xmin, nrow(det)), det$lat), lonlat = lonlat) / cell_size_y

  coa_adj <-
    data.frame(
      x = distances_y,                  #lat
      y = det$dep / zscale,             #depth
      z = abs(distances_x) - ncol(ras)  #lon
    ) %>%
    as.matrix()

  H.pi <- ks::Hpi(coa_adj, binned = TRUE) * mul
  fhat <- ks::kde(coa_adj, H = H.pi)

  ks::plot.kde(
    fhat,
    add = TRUE,
    axes = F,
    box = F,
    lit = F,
    shininess = 128,
    specular = "black",
    fog = F,
    smooth = 2,
    ...
  )
}
