## Plot Bathymetry 


plot_bath<- function (hillshade, heightmap, zscale = 1,
                      water = FALSE, waterdepth = 0, watercolor = "lightblue", 
                      wateralpha = 0.5, theta = 45, phi = 45, fov = 0, zoom = 1, 
                      background = "white", windowsize = c(600, 600), ...) {
  require(rayshader)
  argnameschar = unlist(lapply(as.list(sys.call()), as.character))[-1]
  argnames = as.list(sys.call())
  if (!is.null(attr(heightmap, "rayshader_data"))) {
    if (!("zscale" %in% as.character(names(argnames)))) {
      if (length(argnames) <= 3) {
        zscale = 200
        message("`montereybay` dataset used with no zscale--setting `zscale=200` for a realistic depiction. Lower zscale (i.e. to 50) in `plot_3d` to emphasize vertical features.")
      }
      else {
        if (!is.numeric(argnames[[4]]) || !is.null(names(argnames))) {
          if (names(argnames)[4] != "") {
            zscale = 200
            message("`montereybay` dataset used with no zscale--setting `zscale=200` for a realistic depiction. Lower zscale (i.e. to 50) in `plot_3d` to emphasize vertical features.")
          }
        }
      }
    }
  }
  if (any(hillshade > 1 || hillshade < 0)) {
    stop("Argument `hillshade` must not contain any entries less than 0 or more than 1")
  }
  flipud = function(x) {
    x[nrow(x):1, ]
  }
  if (class(hillshade) == "array") {
    hillshade[, , 1] = flipud(hillshade[, , 1])
    hillshade[, , 2] = flipud(hillshade[, , 2])
    hillshade[, , 3] = flipud(hillshade[, , 3])
  }
  if (class(hillshade) == "matrix") {
    hillshade = hillshade[, ncol(hillshade):1]
  }
  if (is.null(heightmap)) {
    stop("heightmap argument missing--need to input both hillshade and original elevation matrix")
  }
  
  if (water) {
    if (watercolor == "imhof1") {
      watercolor = "#defcf5"
    }
    else if (watercolor == "imhof2") {
      watercolor = "#337c73"
    }
    else if (watercolor == "imhof3") {
      watercolor = "#4e7982"
    }
    else if (watercolor == "imhof4") {
      watercolor = "#638d99"
    }
    else if (watercolor == "desert") {
      watercolor = "#caf0f7"
    }
    else if (watercolor == "bw") {
      watercolor = "#dddddd"
    }
    else if (watercolor == "unicorn") {
      watercolor = "#ff00ff"
    }
  }
  tempmap = tempfile()
  save_png(hillshade, tempmap)
  rgl.surface(1:nrow(heightmap), -(1:ncol(heightmap)), 
              heightmap[, ncol(heightmap):1]/zscale, 
              texture = paste0(tempmap, ".png"), lit = FALSE, ...)
  rgl.bg(color = background)
  rgl.viewpoint(zoom = zoom, phi = phi, theta = theta, fov = fov)
  par3d(windowRect = c(0, 0, windowsize))
  if (water) {
    triangles3d(matrix(c(nrow(heightmap), 1, nrow(heightmap), 
                         waterdepth, waterdepth, waterdepth, -ncol(heightmap), 
                         -1, -1), 3, 3), color = watercolor, alpha = wateralpha, 
                depth_mask = TRUE, front = "fill", back = "fill", depth_test = "lequal")
    triangles3d(matrix(c(1, 1, nrow(heightmap), waterdepth, 
                         waterdepth, waterdepth, -ncol(heightmap), -1, -ncol(heightmap)), 
                       3, 3), color = watercolor, alpha = wateralpha, depth_mask = TRUE, 
                front = "fill", back = "fill", depth_test = "lequal")
  }
}