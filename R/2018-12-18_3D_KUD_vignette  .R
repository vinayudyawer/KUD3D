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

abacusPlot(ATTdata, facet=T)

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

tag <-
  tag_utm %>%
  as_Spatial() %>% as.tibble() %>%
  transmute(X = coords.x1,
            Y = coords.x2,
            depth = Sensor.Value.coa) %>%
  as.matrix

H.pi <- Hpi(tag, binned = TRUE) * 3
fhat <- kde(tag, H = H.pi)

## Create function to calculate 50% and 95% 3DKUD area from fhat
vol3d<-function(fhat, cont=50){
  ct<-contourLevels(fhat, cont=cont, approx=TRUE)
  vol.voxel<- prod(sapply(fhat$eval.points, diff)[1,])
  no.voxel<- sum(fhat$estimate>ct)
  no.voxel*vol.voxel
}

vol3d(fhat, cont = 50) ## in m3
vol3d(fhat, cont = 95)


################################
## Plotting 3D KUD and bathymetry
bath <- crop(bath_ras_utm, extent(stat) + 500) *-1

# ## reconfigure bathymetry data for 3D plotting (** to correct mirrored plotting in rayshader)
bath_mat <-
  as.matrix(bath) %>%
  apply(., 2, rev)

## Set depth exaggeration
depth_exaggeration = 1.5

## plot bathymetry
bath_mat %>%
  sphere_shade(texture = "desert") %>%
  add_shadow(ray_shade(bath_mat, zscale = 1/depth_exaggeration), 0.1) %>%
  add_shadow(ambient_shade(bath_mat, zscale = 1/depth_exaggeration), 0.1) %>%
  plot_3d(
    bath_mat,
    water = T,      ## render water
    zscale = 1/depth_exaggeration,
    waterdepth = 0,
    watercolor = "#88DDFF",
    shadow = FALSE,   ## render shadow on the bottom
    solid = T,    ## render solid base
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
    alpha = 0.3              ## transparency of bathymetry
  )

## Divides COA data into dddn, calculates KUD and then plots it on bathymetry
coa_df %>%
  group_by(dddn) %>%
  transmute(lon = X,
            lat = Y,
            dep = Z) %>%
  do(
    add_kud(
      reef = bath,
      det = .,
      zscale = 1,
      cont = c(95, 50),                        ## you can add multiple contours, but the plot might get a bit cluttered
      alphavec = c(0.1, 0.9),
      drawpoints = F,
      size = 1,
      # col.pt="black", colors=c("red","red")
      col.pt = rainbow(4)[as.numeric(.$dddn[1])],
      colors = rep(rainbow(4)[as.numeric(.$dddn[1])], 2)
    )
  )

## add receiver stations
stat %>%
  data.frame() %>%
  transmute(lon = coords.x1,
            lat = coords.x2,
            dep = -Depth) %>%      ## the depth data for receivers if you have it otherwise you can plot them just under the water surface (-1)
  add_points(reef = bath,
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
add_axes(bath,
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



