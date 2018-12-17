## Estimating and Visualising 3D KUD 
## 
## uses VUE export files with station information metadata to calculate COA
## COAs then used to estimate 3DKUD and visualise as multiple formats

## loading required libraries 
sapply(c("chron","ks","rgl","misc3d","maptools","maps","reshape","PBSmapping","brainR","raster","plyr","lubridate","RColorBrewer"), require, character.only = TRUE)

#######################
## Input files
#######################
## detection data from folder with multiple VUE csv exports
tagdata<-llply(list.files("~/...", full.names = T), read.csv,header=TRUE,.progress="text")
summary(tagdata)

## Station information metadata
ll<-CRS("+proj=longlat +datum=WGS84"); utm<-CRS("+init=epsg:32755")
statinfo<-read.csv("~/../statinfo.csv",header=TRUE,sep=",",strip.white=T) 
stat<-statinfo; coordinates(stat)<-~longitude+latitude; projection(stat)<-ll; stat<-spTransform(stat,utm)

## Bathymetry data
dep<-read.csv("~/.../TSVReef100.csv", sep=",", header=T, strip.white=T)
dep[dep$depth>0,"depth"]<-dep[dep$depth>0,"depth"]/10 ## reduce topography of islands (depth exaggeration)
s <- data.frame(X=dep$X, Y=dep$Y, Z=-dep$depth)
xa<-sort(unique(s$X)); ya<- sort(unique(s$Y)); za<-as.matrix(cast(s,X~Y, value="Z")) 
bath<-rasterFromXYZ(s, crs=utm) ## lower z limit (bathymetry)
ss<-rasterFromXYZ(data.frame(X=s$X, Y=s$Y, Z=0), crs=proj) ## upper z limit (sea surface)

#######################
## COA estimation
#######################
tag<-tagdata[[1]]

## merge tagdata with station information to get latitude and longitude for each detection on array
data<-merge(data.frame(dt=as.POSIXct(strptime(tag[,grep("*Date*",colnames(tag))], format="%Y-%m-%d %H:%M:%S"), tz="GMT"),
                       tagid=sapply(strsplit(as.character(tag$Transmitter),"-"), function(x){as.numeric(x[3])}),
                       serial=sapply(strsplit(as.character(tag$Receiver),"-"), function(x){as.numeric(x[2])}),
                       depth=tag$Sensor.Value), statinfo@data[c("serial","latitude","longitude","stationreef","reefkm","rec_depth")], by="serial", all.x=T)

step<-3600 ## COA timestep in seconds
ex<-seq(from=trunc(min(data$dt, na.rm=TRUE), "day"), to=trunc(max(data$dt, na.rm=TRUE), "day")+86400, by=step)
data$DateTime<-cut(data$dt, breaks=ex) ## add column in dataframe with hourly bins
COA<-ddply(data, "DateTime", summarize, 
             tagid=unique(data$tagid), 
             meanlat=mean(latitude), 
             meanlon=mean(longitude),
             Z=-mean(depth),
             meanreefkm=mean(reefkm))
COA<-COA[!is.na(COA$meanlat),]
cenac<-COA; coordinates(cenac)<-~meanlon+meanlat; projection(cenac)<-ll
pr<-spTransform(cenac, utm); COA[c("X","Y")]<-pr@coords
plot3d(COA$X, COA$Y, COA$Z, aspect=c(1,1,0.3))

#######################
## 3D KUD calculations
#######################
a<-as.matrix(data.frame(X=COA$X, Y=COA$Y, Z=COA$Z))
H.pi <- Hpi(a,binned=TRUE)*3; fhat <- kde(a, H=H.pi)

## Create function to calculate 50% and 95% 3DKUD area from fhat
vol3d<-function(fhat, cont=50){
  ct<-contourLevels(fhat, cont=cont, approx=TRUE)
  vol.voxel<- prod(sapply(fhat$eval.points, diff)[1,])
  no.voxel<- sum(fhat$estimate>ct)
  no.voxel*vol.voxel
}

#######################
## plotting 3DKUD in R
#######################
#### Find x,y,z plotting limits ###
open3d(windowRect=c(46,44,1283,792))
view3d(theta = 0, phi =-45, fov=60, zoom=1)
with(COA, plot3d(X, Y, Z, aspect=c(1,1,0.35), pch=20,cex=0.1, col="blue", xlab="", ylab="", zlab="", axes=F, add=F))
width<-max((par3d("bbox")[2]-par3d("bbox")[1])/2,(par3d("bbox")[4]-par3d("bbox")[3])/2)
xcen<-(par3d("bbox")[2]+par3d("bbox")[1])/2; ycen<-(par3d("bbox")[4]+par3d("bbox")[3])/2
xlim<-c(xcen-(width*2), xcen+(width*2)); ylim<-c(ycen-(width*2), ycen+(width*2)); zlim<-c(par3d("bbox")[5:6])
### subset surfaces to fit plotting limits 
xaz<-xa[xa<xlim[2]&xa>xlim[1]]; yaz<-ya[ya<ylim[2]&ya>ylim[1]]; zaz<-subset(za, xa<xlim[2]&xa>xlim[1], ya<ylim[2]&ya>ylim[1])
### subset recievers to fit plotting limits
recievers<-stat[stat@coords[,1]<xlim[2]&stat@coords[,1]>xlim[1]&stat@coords[,2]<ylim[2]&stat@coords[,2]>ylim[1],]
### colour assignment
zlimz <- range(zaz)
zlenz<- zlimz[2]-zlimz[1]+1
colfunc<-colorRampPalette(c("aquamarine","burlywood3","darkslategray4","dodgerblue3","burlywood3"))
colfuncL<-colorRampPalette(c("black","darkgrey","tan","tan","tan","burlywood1","burlywood2","burlywood3"))
colorlutz <- colfunc(zlenz)
colorlutLz<-rev(colfuncL(zlenz*1.5))
colz <- colorlutz[zaz-zlimz[1]+1]
colLz<- colorlutLz[zaz-zlimz[1]+1]

#### Plot overall detections and KUD
# open3d(windowRect=c(46,44,1283,792))
# view3d(theta = 0, phi =-50, fov=60, zoom=1)
with(COA, plot3d(X, Y, Z, aspect=c(1,1,0.5), type="n",cex=0.1, col="blue", xlab="", ylab="", zlab="", axes=F))
surface3d(xaz,yaz,-zaz,col=colLz,back="filled", lit=F, shininess=128,specular="black" ,fog=F) # all
surface3d(xaz,yaz,-zaz+0.1,col=colz,alpha=0.3, add=T, lit=F) # bath
surf<-matrix(runif(length(yaz)*length(xaz), -0.1,0), length(yaz),length(xaz)) 
surface3d(xaz,yaz,surf,col='#ADD8E6', alpha=0.15,add=TRUE, lit=F)# sea surface
with(COA, points3d(X,Y,Z, pch=20, col="red", alpha=0.2))
plot(fhat,cont=c(50),colors=c("red") , drawpoints=F, add=T, axes=F, box=F, lit=F, shininess=128,specular="black", fog=F, smooth=2, alpha=0.9)
plot(fhat,cont=c(95),colors=c("red") , drawpoints=F, add=T, axes=F, box=F, lit=F, shininess=128,specular="black", fog=F, smooth=2, alpha=0.3)
plot3d(recievers@coords[,1], recievers@coords[,2], -recievers$rec_depth, add=T, type="s", size=1)

######################################################################
##### Saving output in different formats (leave rgl window open) #####
######################################################################
### .GIF animation
movie3d(spin3d(axis=c(0,0,1), rpm=1),dir="~/Desktop/Output", duration=60, fps=10, movie="Spin animation")

### Export rgl window to a web browser
browseURL(paste("file://", writeWebGL(dir=file.path(tempdir(), "webGL"), width=1100), sep=""))

### save as a WebGL()
writeWebGL(dir="Interactive plot", width=1100, snapshot=T)

### Write to .OBJ so can be uploaded to p3d.in server or pdf document
#filename<-paste("~/Desktop/Output/3D KUD/OBJ files/",tagdata$ID[1],".obj", sep="")
#writeOBJ(filename)

### Write to .PLY and .STL format for 3D printing (combines all objects to one single object)
#filename<-paste("~/Desktop/Output/3D KUD/OBJ files/",tagdata$ID[1],".obj", sep="")
#writePLY(filename)
#writeSTL(filename)

###################################################
#### Interactive html plot
###################################################
### Rescale bathymetry
# Since we are working in projected coordinates and values of 
# easting and northing are above 10000, brainR outputs
# will have trouble plotting the full layers.
# So we rescale the bathymetry, detection data, statinfo data and 3DKUD plots
# so they display properly in interactive form

### rescale bathymetry
mul<-10 # depth multiplier for depth exaggeration
base<-surfaceTriangles(x=xaz-xlim[1],y=yaz-ylim[1],-zaz*mul,color="moccasin", fill=T, material="dull", alpha=0.6, smooth=20)
surf<-matrix(runif(length(yaz)*length(xaz), 0,0),length(yaz),length(xaz)) 
top<- surfaceTriangles(x=xaz-xlim[1],y=yaz-ylim[1],surf,color='#ADD8E6', fill=T, alpha=0.4, material="dull", smooth=2)

### rescale detection data and recalculate 3DKUDs for interactive plot  
COAre<-data.frame(X=(COA$X-xlim[1]),Y=(COA$Y-ylim[1]), Z=COA$Z*mul)
H.rescale <- Hpi(COAre,binned=TRUE)*mul
fhat.rescale<- kde(COAre, H=H.rescale)

### Preview bathymetry and 3dKUD
open3d(windowRect=c(46,44,1283,792))
bg3d("#d9d9d9")
drawScene.rgl(base, lighting=perspLighting)
drawScene.rgl(top,lighting=perspLighting, add=T)
with(COAre, plot3d(X, Y, Z, col="red", xlab="", ylab="", zlab="", axes=F, add=T))
plot(fhat.rescale,cont=c(50,95),colors=c(2,2) , drawpoints=F, add=T, axes=F, box=F, lit=F, shininess=128,specular="black", fog=F, smooth=20, alpha=c(0.4,0.8))

### 3DKUD outputs from rgl package have text embedded in .str or .obj exports
### so we have to create dummy files to encorporate 3DKUD voxels in final interactive model
### Create dummy files
stat<- surfaceTriangles(x=xaz-xlim[1],y=yaz-ylim[1],surf,color=1, fill=T, material="dull")
pt<- surfaceTriangles(x=xaz-xlim[1],y=yaz-ylim[1],surf,color=2, fill=T, material="dull")
kud50<- surfaceTriangles(x=xaz-xlim[1],y=yaz-ylim[1],surf,color=2, fill=T, alpha=0.7, material="dull")
kud95<- surfaceTriangles(x=xaz-xlim[1],y=yaz-ylim[1],surf,color=2, fill=T, alpha=0.2, material="dull")

### setting scene values for 4D plot
scene<- list(base,top,stat,pt,kud50,kud95)

setwd("~/Desktop/Example") ## create empty folder named "Example" in desktop to add interactive model
fnames <- c("base.stl", "top.stl","stat.stl",
            "pt.stl","kud50.stl","kud95.stl")
outfile <-  "index.html"
write4D(scene=scene, fnames=fnames, outfile=outfile, visible=c(T,T,T,F,F,F),
        standalone=TRUE, rescale=TRUE, writefiles=T, captions=c("Bathymetry", "Surface", "Recievers", 
                                                                "Center of Actity positions","Core activity space (50%3DKUD)","Extent of activity space (95%3DKUD)"))
rgl.close()

### create reciever, 3DKUD .stl files to replace dummy files in 'Example' folder
### you may get some warning messages after executing writeSTL() function, but you can ignore these,
### R is just warns that the text associated with rgl outputs will not be included in the .stl outputs
setwd("~/Desktop/Example")
plot3d(recievers@coords[,1]-xlim[1], recievers@coords[,2]-ylim[1], -recievers$rec_depth*mul, type="s", col=1, add=T, size=1)
writeSTL("stat.stl"); rgl.close()

pts<-COAre[!duplicated(COAre),] ## remove duplicates to reduce .stl file size
plot3d(pts$X, pts$Y, pts$Z, type="p", size=0.2, xlab="", ylab="", zlab="", axes=F)
writeSTL("pt.stl"); rgl.close()

plot(fhat.rescale, cont=50, colors=2, approx.cont=T, drawpoints=F, add=F, axes=F, box=F, lit=F, shininess=128,specular="black", fog=F, smooth=20)
writeSTL("kud50.stl", pointRadius=0); rgl.close()
plot(fhat.rescale, cont=95, colors=2, approx.cont=T, drawpoints=F, add=F, axes=F, box=F, lit=F, shininess=128,specular="black", fog=F, smooth=20)
writeSTL("kud95.stl", pointRadius=0); rgl.close()


### Now open the index.html file with a web browser (Chrome or Safari recommended) in the 'Example' folder on your desktop
### you can modify the appearance and background color of final interactive models by modifying .html code in the index.html file
### let me know if you need help with modifying the .html code

