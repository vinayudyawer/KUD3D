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

#######################
## Input files
#######################
ll<-CRS("+proj=longlat +datum=WGS84"); utm<-CRS("+init=epsg:32755")

## Station information
load("data/statinfo.RData")

## Tag metadata information
load("data/taginfo.RData")

## Detection data
load("data/tagdata.RData")

## Calculate COA using the Animal Tracking Toolbox
ATTdata<-setupData(Tag.Detections = tagdata,
                   Tag.Metadata = taginfo,
                   Station.Information = statinfo,
                   source = "VEMCO")

# abacusPlot(ATTdata, facet=T)

COAdata <- COA(ATTdata = ATTdata,
               timestep = 60)

## Bathymetry data
load("data/bath.RData")
gbr30<-rasterFromXYZ(bath, crs=ll)

#######################
## Project input data
#######################

tag_ll <-
  COAdata %>%
  filter(Tag.ID %in% 13792) %>%
  st_as_sf(coords = c("Longitude.coa", "Latitude.coa"),
           crs = 4326)

tag_utm <- st_transform(tag_ll, crs = 32755)

lodestone_ll <- crop(gbr30,
                     statinfo %>%
                       filter(installation_name %in% "Lodestone") %>%
                       st_as_sf(coords=c("station_longitude", "station_latitude"), crs=4326) %>%
                       as_Spatial() %>% extent() + 0.05)

lodestone_utm <-
  gbr30 %>%
  projectRaster(., crs=utm) %>%
  crop(.,
       statinfo %>%
         filter(installation_name %in% "Lodestone") %>%
         st_as_sf(coords=c("station_longitude", "station_latitude"), crs=4326) %>%
         st_transform(crs = 32755) %>%
         as_Spatial() %>% extent() + 5000)

#######################
## 3D KUD calculations
#######################

### Process tag data to
# 1) adjust depth measures based on bathymetry
# 2) calculate spline trajectories

tag <-
  tag_utm %>%
  as_Spatial() %>%
  as.tibble() %>%
  transmute(X = coords.x1,
            Y = coords.x2,
            depth = -Sensor.Value.coa,
            dt = TimeStep.coa) %>%
  mutate(upper = 0,
         lower =
           raster::extract(
           x = lodestone_utm,
           y = as_Spatial(tag_utm),
           method='bilinear'),
         Z = case_when(
           depth <= lower ~  lower-(depth-lower),
           depth > lower ~ depth
         ),
         Zadj = ((lower-upper) * (Z-min(Z))/(max(Z)-min(Z))) + upper
         )

plot3d(x=tag$X, y=tag$Y, z=tag$depth, col=1)
points3d(x=tag$X, y=tag$Y, z=tag$Z, col=2)
points3d(x=tag$X, y=tag$Y, z=tag$Zadj, col=3)


COA[COA$Z>COA$bathz,"Z"]<-COA[COA$Z>COA$bathz,"bathz"]-(COA[COA$Z>COA$bathz,"Z"]-COA[COA$Z>COA$bathz,"bathz"])
COA$Zadj<-with(COA, ((bathz-ssz)*(Z-min(Z))/(max(Z)-min(Z)))+ssz)



H.pi <- Hpi(tag, binned = TRUE) * 3
fhat <- kde(tag, H = H.pi)

source("R/vol3d.R")
vol3d(fhat, cont = 50) ## in m3
vol3d(fhat, cont = 95)

################################
## Plotting 3D KUD and bathymetry

# ## reconfigure bathymetry data for 3D plotting (** to correct mirrored plotting in rayshader)
bath_mat <-
  as.matrix(lodestone_utm) %>%
  apply(., 2, rev)

## Set depth exaggeration
depth_exaggeration = 1.2

## plot bathymetry
bath_mat %>%
  sphere_shade(texture = "imhof1") %>%
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

## Plotting using our modified plot_bath() function to control transparency
source("R/plot_bath.R")
bath_mat %>%
  sphere_shade(texture = "imhof1") %>%
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
    alpha = 0.3              ## transparency of bathymetry
  )

## Divides COA data into dddn, calculates KUD and then plots it on bathymetry
source("R/add_fkud.R")
tag %>%
  as_tibble() %>%
  transmute(lon = X,
            lat = Y,
            dep = -depth) %>%
  add_fkud(
    ras = lodestone_utm,
    det = .,
    zscale = 1,
    cont = c(95, 50),                     ## you can add multiple contours, but the plot might get a bit cluttered
    alphavec = c(0.1, 0.9),
    drawpoints = T,
    size = 1,
    col.pt="black", colors=c("red","red")
  )

## add receiver stations
source("R/add_points.R")
statinfo %>%
  filter(installation_name %in% "Lodestone") %>%
  st_as_sf(coords=c("station_longitude", "station_latitude"), crs=4326) %>%
  st_transform(crs=32755) %>%
  as_Spatial() %>%
  data.frame() %>%
  transmute(lon = coords.x1,
            lat = coords.x2,
            dep = -10) %>%      ## the depth data for receivers if you have it otherwise you can plot them just under the water surface (-1)
  add_points(ras = lodestone_utm,
             det = .,
             zscale = 1,
             size = 6,
             col = "black")

## add legend
legend3d("topright",
         legend = c(levels(coa_df$dddn), NA , "Receiver stations"),
         border = c(rainbow(4), NA, NA),
         fill = c(rainbow(4), NA, NA),
         pch = c(NA,NA,NA,NA,NA,19),
         col = c(NA,NA,NA,NA,NA,1),
         cex=1, inset=c(0.06))

## add axes
source("R/add_axes.R")
add_axes(lodestone_utm,
         zscale = 1/depth_exaggeration,
         axis.col = grey(0.5))

######################################################################
##### Saving output in different formats (leave rgl window open) #####
######################################################################
### Save output as a .png
snapshot3d(paste0("3Dkud_",data$tag[1],".png"))

### .GIF animation
movie3d(spin3d(axis=c(0,1,0), rpm=1),dir="~/Desktop/Output", duration=60, fps=10, movie="Spin animation")

### Export rgl window to a web browser
browseURL(paste("file://", writeWebGL(dir=file.path(tempdir(), "webGL")), sep=""))

### save as a WebGL()
# writeWebGL(dir="Interactive plot", snapshot=T)

### Write to .OBJ so can be uploaded to p3d.in server or pdf document
#filename<-paste("~/Desktop/Output/3D KUD/OBJ files/",tagdata$ID[1],".obj", sep="")
#writeOBJ(filename)

### Write to .PLY and .STL format for 3D printing (combines all objects to one single object)
#filename<-paste("~/Desktop/Output/3D KUD/OBJ files/",tagdata$ID[1],".obj", sep="")
#writePLY(filename)
#writeSTL(filename)



