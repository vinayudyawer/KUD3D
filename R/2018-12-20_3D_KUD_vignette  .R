## Plotting 3DKUD Heron and Opal Reef sites

## Rayshader: package to help plot 3d bathymetry data
# devtools::install_github("tylermorganwall/rayshader")

## loading required libraries
sapply(c("ks",
         "rgl",
         "raster",
         "tidyverse",
         "lubridate",
         "sf",
         "rayshader",
         "VTrack"),
       require, character.only = TRUE)

## CRS
ll<-CRS("+proj=longlat +datum=WGS84")
utm_heron<-CRS("+init=epsg:32755")
utm_ningaloo<-CRS("+init=epsg:28348")

###########################
## GPS/Satellite telemetry
###########################
## Input and process data

## Detection data
load("data/GPSdata.RData")

shark_ll <-
  GPSdata %>%
  st_as_sf(coords=c("Longitude","Latitude"), crs=4326)

shark_utm <-
  shark_ll %>%
  st_transform(crs=28348)

## Bathymetry data for Ningaloo Reef
load("data/ningaloo_bath.RData")

ningaloo_ll<- rasterFromXYZ(ningaloo_bath, crs=ll)
ningaloo_utm <-
  projectRaster(ningaloo_ll, crs=utm_ningaloo) %>%
  crop(.,
       shark_utm %>%
         as_Spatial() %>%
         extent() + 500)

#######################
## 3D KUD calculations
#######################
kud_df<-
  shark_utm %>%
  as_Spatial() %>%
  as_tibble() %>%
  transmute(X = coords.x1,
            Y = coords.x2,
            Z = - Depth,
            dt = Date.Time)

H.pi <- Hpi(as.matrix(kud_df[1:3]), binned = TRUE)
fhat <- kde(as.matrix(kud_df[1:3]), H = H.pi)

source("R/vol3d.R")
vol3d(fhat, cont = 50) ## in m3
vol3d(fhat, cont = 95)

################################
## Plotting 3D KUD and bathymetry
## Set depth exaggeration
depth_exaggeration <- 0.1

## reconfigure bathymetry data for 3D plotting (** to correct mirrored plotting in rayshader)
bath_mat <-
  as.matrix(ningaloo_utm) %>%
  apply(., 2, rev)


bath_mat %>%
  sphere_shade(texture = "desert") %>%
  add_shadow(ray_shade(bath_mat, zscale = 1/depth_exaggeration), 0.1) %>%
  add_shadow(ambient_shade(bath_mat, zscale = 1/depth_exaggeration), 0.1) %>%
  plot_3d(
    bath_mat,
    baseshape = "rectangle",
    water = T,      ## render water
    zscale = 1/depth_exaggeration,
    wateralpha = 0.2,
    waterlinecolor = "white",
    waterlinealpha = 0.5,
    windowsize = c(1200, 700),  ## Size of window
    theta = 80,               ## Play around with the theta, phi, fov and zoom to orient the plot (or you can adjust it manually)
    phi = 20,
    fov = 60,
    zoom = 0.8
  )

## Plot 3DKUD ontop of bathymetry
source("R/add_fkud.R")
kud_df %>%
  transmute(lat = Y,
            lon = X,
            dep = Z) %>%
  add_fkud(
    ras = ningaloo_utm,
    det = .,
    zscale = 1 / depth_exaggeration,
    cont = c(95, 50),
    alphavec = c(0.1, 0.9),
    drawpoints = T,
    size = 1,
    col.pt = "black",
    colors = c("red","red")
  )

## add axes
source("R/add_axes.R")
add_axes(ningaloo_utm,
         zscale = 1/depth_exaggeration,
         axis.col = grey(0.5))



##########################
## Passive telemetry data
##########################
## Input and process files

## Station information
load("data/statinfo.RData")

stat_ll <-
  statinfo %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

stat_utm <-
  stat_ll %>%
  st_transform(crs = 32755)

## Centre of activity data
load("data/COAdata.RData")

trout_ll <-
  COAdata %>%
  st_as_sf(coords = c("Longitude.coa", "Latitude.coa"), crs = 4326)

trout_utm <-
  trout_ll %>%
  st_transform(crs = 32755)

## Bathymetry data for Heron Island Reef
load("data/heron_bath.RData")

heron_ll<-
  rasterFromXYZ(heron_bath, crs=ll) %>%
  raster::disaggregate(x=., fact=3, method='bilinear')

heron_utm<-
  projectRaster(heron_ll, crs=utm_heron) %>%
  raster::crop(.,
       stat_utm %>%
         as_Spatial() %>%
         extent() + 300)


#######################
## 3D KUD calculations
#######################
kud_df<-
  trout_utm %>%
  as_Spatial() %>%
  as_tibble() %>%
  transmute(X = coords.x1,
            Y = coords.x2,
            Z = - Sensor.Value.coa,
            subset =
              factor(
                case_when(
                  lubridate::hour(TimeStep.coa) %in% c(7:17) ~ "Day",
                  lubridate::hour(TimeStep.coa) %in% c(0:6, 18:23) ~ "Night"
                ), levels=c("Night", "Day")))


H.pi <- Hpi(as.matrix(kud_df[1:3]), binned = TRUE)
fhat <- kde(as.matrix(kud_df[1:3]), H = H.pi)

source("R/vol3d.R")
vol3d(fhat, cont = 50) ## in m3
vol3d(fhat, cont = 95)

################################
## Plotting 3D KUD and bathymetry

## reconfigure bathymetry data for 3D plotting (** to correct mirrored plotting in rayshader)
bath_mat <-
  as.matrix(heron_utm) %>%
  apply(., 2, rev)

## Set depth exaggeration
depth_exaggeration <- 1.5

## Plotting using our modified plot_bath() function to control transparency
source("R/plot_bath.R")
bath_mat %>%
  sphere_shade(texture = "desert") %>%
  add_shadow(ray_shade(bath_mat, zscale = 1/depth_exaggeration), 0.1) %>%
  add_shadow(ambient_shade(bath_mat, zscale = 1/depth_exaggeration), 0.1) %>%
  plot_bath(
    bath_mat,
    water = TRUE,      ## render water surface
    zscale = 1/depth_exaggeration,
    waterdepth = 0,
    watercolor = "#88DDFF",
    wateralpha = 0.2,
    windowsize = c(1200, 700),  ## Size of window
    theta = 80,               ## Play around with the theta, phi, fov and zoom to orient the plot (or you can adjust it manually)
    phi = 20,
    fov = 60,
    zoom = 0.8,
    alpha = 0.4              ## transparency of bathymetry
  )

## Divides COA data into subsets, calculates KUD and then plots it on bathymetry
kud_df %>%
  group_by(subset) %>%
  transmute(lon = X,
            lat = Y,
            dep = Z) %>%
  do(
    add_fkud(
      ras = heron_utm,
      det = .,
      zscale = 1,
      cont = c(95, 50),                        ## you can add multiple contours, but the plot might get a bit cluttered
      alphavec = c(0.1, 0.7),
      drawpoints = F,
      size = 1,
      col.pt = rainbow(2)[as.numeric(.$subset[1])],
      colors = rep(rainbow(2)[as.numeric(.$subset[1])], 2)
    )
  )

## add receiver stations
source("R/add_points.R")
stat_utm %>%
  as_Spatial() %>%
  data.frame() %>%
  transmute(lon = coords.x1,
            lat = coords.x2,
            dep = - Depth) %>%      ## the depth data for receivers if you have it otherwise you can plot them just under the water surface (-1)
  add_points(ras = heron_utm,
             det = .,
             zscale = 1,
             size = 6,
             col = "black")

## add legend
legend3d("topright",
         legend = c(levels(kud_df$subset), NA , "Receiver stations"),
         border = c(rainbow(2), NA, NA),
         fill = c(rainbow(2), NA, NA),
         pch = c(NA,NA,NA,19),
         col = c(NA,NA,NA,1),
         cex=1, inset=c(0.06))

## add axes
add_axes(heron_utm,
         zscale = 1/depth_exaggeration,
         axis.col = grey(0.5))


######################################################################
##### Saving output in different formats (leave rgl window open) #####
######################################################################
### Save output as a .png
snapshot3d("RGLsnapshot.png")

### .GIF animation
movie3d(
  spin3d(axis = c(0, 1, 0), rpm = 1),
  dir = "~/Desktop/Output",
  duration = 60,
  fps = 10,
  movie = "Spin animation"
)

### Export rgl window to a web browser
browseURL(
  paste("file://", writeWebGL(dir = file.path(tempdir(), "webGL")), sep = ""))

### save as a WebGL()
# writeWebGL(dir="Interactive plot", snapshot=T)

### Write to .OBJ so can be uploaded to p3d.in server or pdf document
#filename<-paste("~/Desktop/Output/3D KUD/OBJ files/",tagdata$ID[1],".obj", sep="")
#writeOBJ(filename)

### Write to .PLY and .STL format for 3D printing (combines all objects to one single object)
#filename<-paste("~/Desktop/Output/3D KUD/OBJ files/",tagdata$ID[1],".obj", sep="")
#writePLY(filename)
#writeSTL(filename)



