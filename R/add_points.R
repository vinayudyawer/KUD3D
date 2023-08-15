#' Add points or lines to a 3D rayshader plot
#'
#' @description Helper function to add points or lines a 3D rayshader plot
#'
#' @param ras raster object with bathymetry. This is the raster used to create the hillshade nad heightmap of the rayshader 3D plot
#'   tag detection data, metadata and station information
#' @param det ...
#' @param zscale scaling factor for the depth, indicating depth exaggeration
#' @param lonlat ...
#' @param col colour of point or line
#' @param alpha transparency of point or line
#' @param line logical if to plot lines instead of points
#' @param ... additional arguments for \code{\link{points3d}} or \code{\link{lines3d}} function depending on what your plotting
#'
#' @return Adds points or lines to an existing 3D rayshader plot
#'
#' @import raster
#' @import rgl
#' @importFrom raster extent
#' @importFrom raster pointDistance
#' @importFrom rgl points3d
#' @importFrom rgl lines3d
#'
#' @export add_points
#'
#' @examples
#' TBD
#'

add_points <-
  function(ras,
           det,
           zscale,
           lonlat = FALSE,
           col = "red",
           alpha = 0.8,
           line = F,
           ...) {
    e <- extent(ras)
    cell_size_x <-
      raster::pointDistance(c(e@xmin, e@ymin), c(e@xmax, e@ymin), lonlat = lonlat) / ncol(ras)
    cell_size_y <-
      raster::pointDistance(c(e@xmin, e@ymin), c(e@xmin, e@ymax), lonlat = lonlat) / nrow(ras)
    distances_x <-
      raster::pointDistance(c(e@xmin, e@ymin), cbind(det$lon, rep(e@ymin, nrow(det))), lonlat = lonlat) / cell_size_x
    distances_y <-
      raster::pointDistance(c(e@xmin, e@ymin), cbind(rep(e@xmin, nrow(det)), det$lat), lonlat = lonlat) / cell_size_y

    rgl::points3d(
      x = distances_y - (nrow(ras)/2),                       #lat
      y = det$dep / zscale,                  #depth
      z = abs(distances_x) - (ncol(ras)/2),      #lon
      color = col,
      alpha = alpha,
      ...
    )
    if (line) {
      rgl::lines3d(
        x = distances_y - (nrow(ras)/2),                     #lat
        y = det$dep / zscale,                #depth
        z = abs(distances_x) - (ncol(ras)/2),    #lon
        color = col,
        alpha = alpha,
        ...
      )
    }

  }
