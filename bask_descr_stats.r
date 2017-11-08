# some bask descriptive stats

meta <- read.table('~/ebs/Data/BaskingSharks/batch/bask_metadata_v2.csv', header=T, sep=',', blank.lines.skip = F)
load('~/ebs/Data/BaskingSharks/batch/bask_results_v4.rda')
library(raster); library(fields); library(rgdal); library(trip)

meta$dist <- NA
meta$minlat <- NA
df <- bask.res[[3]]
u.ptts <- unique(df$ptt)

for (i in 1:length(u.ptts)){
  locs.i <- df[which(df$ptt == u.ptts[i]),]
  locs.i$date <- as.POSIXct(locs.i$date)
  # convert track to SpatialPointsDataFrame
  coordinates(locs.i) <- ~lon + lat
  proj4string(locs.i)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  
  # then to trip object
  tr <- trip(locs.i, c('date','ptt'))
  meta$dist[which(meta$PTT %in% u.ptts[i])] <- sum(trackDistance(tr))
  meta$minlat[which(meta$PTT %in% u.ptts[i])] <- min(locs.i$lat, na.rm=T)
  
}

shorts <- meta$PTT[which(meta$DeployDuration < 300)]
shorts <- shorts[which(shorts %in% u.ptts)]

for(i in 1:length(shorts)){
  plot(df$lon[which(df$ptt %in% shorts[i])], df$lat[which(df$ptt %in% shorts[i])], type='l')
  world(add=T, fill=T); title(paste(shorts[i],' ', meta$DeployDuration[which(meta$PTT %in% shorts[i])], ' ', meta$TL_cm[which(meta$PTT %in% shorts[i])]))
  invisible(readline(prompt="Press [enter] to continue"))
  
}

zvec <- c()
for (i in 1:length(profs)){
  if(lubridate::month(profs[[i]]$Date) > 8 | lubridate::month(profs[[i]]$Date) < 6){
    zvec <- c(zvec, min(profs[[i]]$Depth, na.rm=T))
  }
}

meta$surf_freq <- NA
for (ii in 1:nrow(meta)){
  idx <- which(unlist(lapply(profs, FUN=function(x) x$Ptt[1] == meta$PTT[ii])))
  s <- c()
  for (t in idx){
    s <- c(s, min(profs[[t]]$Depth, na.rm=T))
  }
  meta$surf_freq[ii] <- length(which(s > 200)) / length(idx)
}

range(meta$surf_freq, na.rm=T)
mean(meta$surf_freq, na.rm=T)
sd(meta$surf_freq, na.rm=T)
