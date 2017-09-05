# FIGURE 1 - BASKING SHARK MS
load('~/ebs/Data/BaskingSharks/batch/bask_results.rda')
library(raster); library(fields); library(rgdal)

# bask.res is a list of length 3
# [[1]] is named list (all ptts) with allRD, behavRD, and df of track for each
# [[2]] is raster stack of allRD for each ptt (each has own layer)
# [[3]] is df of all tracks together

r.pts <- rasterToPoints(bask.res[[2]], spatial=TRUE)

lon <- r.pts@coords[,1]
lat <- r.pts@coords[,2]
rm(r.pts)
grid <- HMMoce::meshgrid(lon, lat)

# calc semi minor axis based on longitude error
slon.sd <- locs$Error.Semi.minor.axis / 1000 / 111 #semi minor axis
L.light.lon <- stats::dnorm(t(g1$X), locs$Longitude, slon.sd) # Longitude data
slat.sd <- locs$Error.Semi.major.axis / 1000 / 111 #semi major axis
L.light.lat <- stats::dnorm(t(g1$Y), locs$Latitude, slat.sd)

#image.plot(g$lon[1,],g$lat[,1],L.light.lat*L.light.lon)

# create the ellipse by multiplying lat * lon error
L <- raster::flip(raster::raster(t(L.light.lat * L.light.lon), xmn = min(lon1), 
                                 xmx = max(lon1), ymn = min(lat1), ymx = max(lat1)), direction = 'y')



