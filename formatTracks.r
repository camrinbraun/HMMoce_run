# Compare "Known" Positions to Other Traditional Outputs


formatTracks <- function(trackDir, ptt, gpeNo, knock=FALSE){

  # lazy load workspace file to get objects that aren't written indiv to disk
  # convert .RData -> .rdb/.rdx
  print('If you receive a database corrupt error, restart your R session. Try .rs.restartR()')
  e <- local({load(paste(trackDir, '/', ptt, '_traditional_geo.RData', sep='')); environment()})
  tools:::makeLazyLoadDB(e, "New")
  lazyLoad("New")

  strt <- as.Date(paste(sp[3], '-', sp[4], '-', sp[5], sep=''))
  end <- as.Date(paste(ep[3], '-', ep[4], '-', ep[5], sep=''))
  dates <- seq(strt, end, by='day')

  # spot has already been dealt with. just part of the list now.
  spot <- read.table(paste(ptt, '_crawl_track.csv', sep=''), sep=',', header=T)
  spot.dates <- as.Date(spot$date)
  spot <- spot[which(spot.dates <= end),]

  ## then re-format each input track type, as necessary
  # TI: trim to 1 pos daily, lon - 360
  ti <- fit
  ti$most.prob.track[,1] <- ti$most.prob.track[,1] - 360
  tidates <- as.Date(ti$date, origin = '1960-01-01')
  idx <- which(!duplicated(tidates) & tidates <= end & tidates >= strt)
  tidates <- tidates[idx]
  ticoords <- ti$most.prob.track[idx,]
  ti.var <- ti$var.most.prob.track[idx,]

  # TI + B: trim to 1 pos daily
  tib <- read.table(paste(ptt, '_fit_btrack.csv', sep=''), sep=',', header=T)
  if(any(tib$Lon_E > 180)) tib$Lon_E <- tib$Lon_E - 360
  tib <- tib[idx,]

  # KF: trim to 1 pos daily, lon - 360
  kf <- kfit
  kf$most.prob.track[,1] <- kf$most.prob.track[,1] - 360
  kfdates <- as.Date(paste(kf$date[,1], '-', kf$date[,2], '-', kf$date[,3], sep=''))
  kf.idx <- which(!duplicated(kfdates) & kfdates <= end & kfdates >= strt)
  kfdates <- kfdates[kf.idx]
  kfcoords <- kf$most.prob.track[kf.idx,]
  kf.var <- kf$var.most.prob.track[kf.idx,]

  # KF + B:
  kfb <- kfbtrack
  kfb <- kfb[kf.idx,]

  # GPE3
  gpe <- read.table(paste(ptt, '-', gpeNo, '-GPE3.csv', sep=''), sep=',', header = T, skip = 5)
  gpedates <- as.Date(gpe$Date, format='%d-%b-%Y %H:%M:%S')
  gpe.idx <- which(!duplicated(gpedates) & gpedates <= end & gpedates >= strt)
  gpedates <- gpedates[gpe.idx]
  crd.idx <- c(grep('lon', names(gpe), ignore.case=T), grep('lat', names(gpe), ignore.case=T))
  gpecoords <- gpe[gpe.idx, crd.idx]

  # HMMoce
  if(knock){
    hmm <- read.table(paste(ptt, '_HMM_track_knock.csv', sep=''), sep=',', header = T)
  } else{
    hmm <- read.table(paste(ptt, '_HMM_track.csv', sep=''), sep=',', header = T)
  }
  hmm.dates <- as.Date(hmm$date)

  # having a df with one row for each day and cols for each track type would be nice.
  # then would just have to save the index of which coords we're keeping so we can
  # pull corresponding var for CI
  df <- as.data.frame(dates)
  df[,2:19] <- NA
  #print(str(df))
  #print(str(which(dates %in% tidates)))
  df[which(dates %in% tidates), 2:3] <- ticoords
  df[which(dates %in% tidates), 4:5] <- ti.var
  df[which(dates %in% tidates), 6:7] <- tib[,c(8:9)]
  df[which(dates %in% kfdates), 8:9] <- kfcoords
  df[which(dates %in% kfdates), 10:11] <- kf.var
  df[which(dates %in% kfdates), 12:13] <- kfb[,c(8:9)]
  df[which(dates %in% spot.dates), 14:15] <- spot[,c(2:3)]
  df[which(dates %in% gpedates), 16:17] <- gpecoords
  df[which(dates %in% hmm.dates), 18:19] <- hmm[,c(2:3)]

  names(df) <- list('date','ti.lon','ti.lat','ti.var.lon','ti.var.lat','tib.lon','tib.lat',
                    'kf.lon','kf.lat','kf.var.lon','kf.var.lat','kfb.lon','kfb.lat',
                    'spot.lon','spot.lat','gpe.lon','gpe.lat','hmm.lon','hmm.lat')

  df

}
