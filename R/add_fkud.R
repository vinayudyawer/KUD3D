## add fixed KUD to rayshader

add_fkud <- function(reef, det, zscale, lonlat = FALSE, mul=1, ...){
  e <- raster::extent(reef)
  cell_size_x <- raster::pointDistance(c(e@xmin, e@ymin), c(e@xmax, e@ymin), lonlat = lonlat) / ncol(reef)
  cell_size_y <- raster::pointDistance(c(e@xmin, e@ymin), c(e@xmin, e@ymax), lonlat = lonlat) / nrow(reef)
  distances_x <- raster::pointDistance(c(e@xmin, e@ymin), cbind(det$lon, rep(e@ymin, nrow(det))), lonlat = lonlat) / cell_size_x
  distances_y <- raster::pointDistance(c(e@xmin, e@ymin), cbind(rep(e@xmin, nrow(det)), det$lat), lonlat = lonlat) / cell_size_y
  
  coa_adj<-
    data.frame(
      x = distances_y,  #lat
      y = det$dep / zscale,  #depth
      z = abs(distances_x) - ncol(reef)  #lon
    ) %>%
    as.matrix()
  
  H.pi <- Hpi(coa_adj,binned=TRUE) * mul
  fhat <- kde(coa_adj, H=H.pi)
  
  plot(fhat, 
       add = TRUE,
       axes = F, 
       box = F,
       lit = F,
       shininess = 128,
       specular = "black",
       fog = F,
       smooth = 2,
       ...)
}