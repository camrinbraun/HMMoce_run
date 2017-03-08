# Calculate speed histogram from format.locs output track

speedHist <- function(track, plot = TRUE){

  if(length(unique(track$Ptt)) > 1) stop('Error: Multiple PTTs in input track file.')

  track <- track[order(track$Date),]

  # convert track to SpatialPointsDataFrame
  coordinates(track) <- ~lon + lat
  proj4string(track)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

  # then to trip object
  tr <- trip(track, c('Date','Ptt'))

  tdist <- trip::trackDistance(tr, longlat = T)
  timediff <- diff(tr@data$Date) # in seconds

  if(units(timediff) == 'hours') timediff <- timediff * 3600
  if(units(timediff) == 'mins') timediff <- timediff * 60
  #if(units(timediff) == 'seconds') timediff <- unclass(timediff)

  timediff <- unclass(timediff)

  spd <- tdist[2:length(tdist)] / timediff

  h <- hist(spd*1000, xlab = 'Speed (m/s)', main='', xlim=c(0,10), breaks=seq(0,ceiling(max(spd)*1000), by=.25))#, ylim=c(0,50))
  #hist(spd*1000, xlab = 'Speed (m/s)', main='')#, xlim=c(0,30), breaks=seq(0,ceiling(max(spd)), by=1), ylim=c(0,50))
  abline(v=mean(h[[3]]), col='red')

}
