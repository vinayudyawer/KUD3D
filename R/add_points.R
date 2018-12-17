## Add points to rayshader


add_points <- function(reef, det, zscale, lonlat = FALSE, col = "red", alpha = 0.8, line=F, ...){
  e <- raster::extent(reef)
  cell_size_x <- raster::pointDistance(c(e@xmin, e@ymin), c(e@xmax, e@ymin), lonlat = lonlat) / ncol(reef)
  cell_size_y <- raster::pointDistance(c(e@xmin, e@ymin), c(e@xmin, e@ymax), lonlat = lonlat) / nrow(reef)
  distances_x <- raster::pointDistance(c(e@xmin, e@ymin), cbind(det$lon, rep(e@ymin, nrow(det))), lonlat = lonlat) / cell_size_x
  distances_y <- raster::pointDistance(c(e@xmin, e@ymin), cbind(rep(e@xmin, nrow(det)), det$lat), lonlat = lonlat) / cell_size_y
  
  rgl::points3d(
    x = distances_y,  #lat
    y = det$dep / zscale,  #depth
    z = abs(distances_x) - ncol(reef),  #lon
    color = col,
    alpha = alpha,
    ...
  )
  if(line){
    rgl::lines3d(
      x = distances_y,  #lat
      y = det$dep / zscale,  #depth
      z = abs(distances_x) - ncol(reef),  #lon
      color = col,
      alpha = alpha,
      ...
    )
  }
  
}