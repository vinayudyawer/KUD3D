## Add plot axes to rayshader

add_axes <- function(reef, zscale, axis.col = 1, ...){
  
  ## define axis labs (ll and utm)
  llabx <-
    pretty(range(coordinates(projectRaster(
      reef, crs = CRS("+proj=longlat +datum=WGS84")
    ))[, 1]))
  llaby <-
    pretty(range(coordinates(projectRaster(
      reef, crs = CRS("+proj=longlat +datum=WGS84")
    ))[, 2]))
  
  labz <- pretty(c(0, min(values(reef))))
  
  ## estimate positions for labels
  atx <- seq(-ncol(reef), 0, len=length(llabx))
  aty <- seq(0, nrow(reef), len=length(llaby))
  atz <- labz / zscale
  
  ## add axes
  bbox3d(xat = aty, xlab = llaby,
         yat = atz, ylab = labz,
         zat = atx, zlab = llabx,
         col = axis.col, marklen = 30)
}