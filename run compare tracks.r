
#spotList <- c(141270, 141268, 141267, 141261, 121420)
satList <- c(141254, 141256, 141257, 141259, 121325)

source('~/Documents/WHOI/RCode/HMMoce_run/formatTracks.r')
source('~/Documents/WHOI/RCode/HMMoce_run/compareTracks.r')

gpeNo <- c(5,1,1,4,1)
path <- '~/Documents/WHOI/RCode/HMMoce_run/data/'

for (ii in 1:length(satList)){
  setwd(paste(path, satList[ii], '/', sep=''))
  df <- formatTracks(trackDir = paste(path, satList[ii], '/', sep=''), ptt = satList[ii], gpeNo = gpeNo[ii])
  res <- compareTracks(df)
}

save(res, file=paste(satList[ii], '_res.RData'))

#df <- formatTracks(trackDir, ptt, 1)
mean(res[[3]]$hmm, na.rm=T); mean(res[[3]]$gpe, na.rm=T)
res[[1]][6]; res[[1]][5]#rmse lon
res[[2]][6]; res[[1]][5] #rmse lat


for (ii in 4:length(satList)){
  setwd(paste(path, satList[ii], '/', sep=''))
  df <- read.table(paste(satList[ii], '_allTracks.csv', sep=''), sep=',', header=T)
  #df <- formatTracks(trackDir = paste(path, satList[ii], '/', sep=''), ptt = satList[ii], gpeNo = gpeNo[ii])
  latlow <- floor(min(df$spot.lat)) - 5; lathigh <- ceiling(max(df$spot.lat)) + 5
  lonlow <- floor(min(df$spot.lon)) - 5; lonhigh <- ceiling(max(df$spot.lon)) + 5
  plot(df$spot.lon, df$spot.lat, type='l', xlim=c(lonlow, lonhigh), ylim=c(latlow, lathigh))
  world(add=T, fill=T, col='grey60')
  lines(df$gpe.lon, df$gpe.lat, col='blue')
  lines(df$hmm.lon, df$hmm.lat, col='red')
  write.table(df, file=paste(satList[ii], '_allTracks.csv', sep=''), sep=',', col.names=T, row.names=F)

}

#===============
## END

# load SPOT data
argos <- read.table('~/Documents/WHOI/RData/sharkSiteData/AllArgosData.csv', sep=',', header=T)
argos$date <- as.POSIXct(argos$date, format='%Y-%m-%d %H:%M:%S')

lydia <- read.table('~/Documents/WHOI/Data/WhiteSharks/2013/121325/121420-Locations.csv', sep=',', header = T)
lydia <- lydia[,c(2,4,7:8,6)]
names(lydia) <- names(argos)
lydia$date <- as.POSIXct(lydia$date, format='%m/%d/%y %H:%M')
argos <- rbind(argos, lydia)
argos <- argos[which(argos$ptt != 141264),]


for (i in 1:length(satList)){
  argos$ptt[which(argos$ptt == spotList[i])] <- satList[i]
}

for (i in 1:length(satList)){
  argos.i <- argos[which(argos$ptt == satList[i]),]

  # remove duplicate date-times
  argos.i <- argos.i[which(!duplicated(argos.i$date)),]

  # convert track to SpatialPointsDataFrame
  coordinates(argos.i) <- ~lon + lat
  proj4string(argos.i)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

  # then to trip object
  tr <- trip(argos.i, c('date','ptt'))

  # filter at 2 m/s or 7.2 km/hr
  sf2 <- speedfilter(tr, max.speed=7.2)

  # and subset track based on logical output of speed filter
  #tr.f <- tr[sf2,]
  argos.i <- data.frame(tr[sf2,])[,c(1:4)]

  if (i == 1){
    argos.new <- argos.i
  } else{
    argos.new <- rbind(argos.new, argos.i)
  }

}

