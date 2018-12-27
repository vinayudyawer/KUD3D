#' Add axes to 3D rayshader plot
#'
#' @description Helper function to add axes to a 3D rayshader plot
#'
#' @param ras raster object with bathymetry. This is the raster used to create the hillshade nad heightmap of the rayshader 3D plot
#'   tag detection data, metadata and station information
#' @param zscale scaling factor for the depth, indicating depth exaggeration
#' @param axis.col color of axis and axis text, defaulted to 1
#' @param ... additional arguments for \code{\link{bbox3d}} function
#'
#' @return Adds axes to an existing 3D rayshader plot
#'
#' @import raster
#' @import sp
#' @import rgl
#' @importFrom raster coordinates
#' @importFrom raster projectRaster
#' @importFrom raster values
#' @importFrom sp CRS
#' @importFrom rgl bbox3d
#'
#' @export add_axes
#'
#' @examples
#' TBD
#'

add_axes <- function(ras, zscale, axis.col = 1, ...){

  ## define axis labs (ll and utm)
  llabx <-
    pretty(range(raster::coordinates(raster::projectRaster(
      ras, crs = sp::CRS("+proj=longlat +datum=WGS84")
    ))[, 1]))
  llaby <-
    pretty(range(raster::coordinates(raster::projectRaster(
      ras, crs = sp::CRS("+proj=longlat +datum=WGS84")
    ))[, 2]))

  labz <- pretty(c(0, min(raster::values(ras))))

  ## estimate positions for labels
  atx <- seq(-ncol(ras), 0, len=length(llabx))
  aty <- seq(0, nrow(ras), len=length(llaby))
  atz <- labz / zscale

  ## add axes
  rgl::bbox3d(xat = aty, xlab = llaby,
         yat = atz, ylab = labz,
         zat = atx, zlab = llabx,
         col = axis.col, marklen = 30, ...)
}

