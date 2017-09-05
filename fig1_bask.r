# FIGURE 1 - BASKING SHARK MS
load('~/ebs/Data/BaskingSharks/batch/bask_results.rda')
library(raster); library(fields); library(rgdal)

# bask.res is a list of length 3
# [[1]] is named list (all ptts) with allRD, behavRD, and df of track for each
# [[2]] is raster stack of allRD for each ptt (each has own layer)
# [[3]] is df of all tracks together

r.pts <- rasterToPoints(bask.res[[2]], spatial=TRUE)

lon.all <- unique(r.pts@coords[,1])
lat.all <- rev(unique(r.pts@coords[,2]))
rm(r.pts)
g.all <- HMMoce:::meshgrid(lon.all, lat.all)

# get a list of SpatialPolygons for each individual
for (i in 1:length(bask.res[[1]])){
  locs <- bask.res[[1]][[i]]$df
  polyList <- list()
  for (b in 1:nrow(locs)){
    # calc axes of ellipse based on x/ydist from getCtr
    slon.sd <- locs$xdist[b]
    L.light.lon <- stats::dnorm(t(g.all$X), locs$lon[b], slon.sd) # Longitude data
    slat.sd <- locs$ydist[b]
    L.light.lat <- stats::dnorm(t(g.all$Y), locs$lat[b], slat.sd)
    
    # create the ellipse by multiplying lat * lon error
    L <-  raster::flip(raster::raster(t(L.light.lat * L.light.lon), xmn = min(lon.all), 
                         xmx = max(lon.all), ymn = min(lat.all), ymx = max(lat.all)), direction='y')
    L.mat <- t(as.matrix(raster::flip(L, direction='y'))) / cellStats(L, 'max')
    ctr.L <- contourLines(lon.all, lat.all, L.mat)
    #idx <- which.min(lapply(ctr, FUN=function(x) which(round(x$level,1) == round(threshold, 1))) == 1)
    ctr <- data.frame(ctr[[1]])
    sp::coordinates(ctr) <- ~x+y
    l1 <- sp::SpatialLines(list(sp::Lines(sp::Line(sp::coordinates(ctr)), "L1")))
    p1 <- gPolygonize(l1)
    polyList[[b]] <- p1
  }
  
  bask.res[[1]][[i]]$polyList <- polyList
  
}



