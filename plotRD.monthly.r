# distr=res$s; track=res$tr; g=res$g; dateVec=res$dateVec

plotRD.monthly <- function(distr, track, ptt, g, dateVec){
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  rd.cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdYlGn')))
  
  p.1 <- apply(distr[1,,,], 1, sum)
  p.2 <- apply(distr[2,,,], 1, sum)
  sv <- -(apply(distr[1,,,], 1, sum) > apply(distr[2,,,], 1, sum)) + 1
  #sv[sv == 0] <- NA
  
  # normalize for each behav at each time point
  norm <- array(NA, dim=dim(distr)[c(3,4,2,1)])
  for (i in 1:dim(norm)[3]){
    norm[,,i,1] <- (distr[1,i,,] / max(distr[1,i,,], na.rm=T)) * p.1[i]
    norm[,,i,2] <- (distr[2,i,,] / max(distr[2,i,,], na.rm=T)) * p.2[i]
  }
  
  all <- apply(norm, 1:3, FUN=function(x) sum(x, na.rm=T))
  
  for (i in 1:dim(all)[3]){
    all[,,i] <- all[,,i] / max(all[,,i], na.rm=T)
  }
  
  track$season <- NA
  
  for (ii in 2004:2012){
    # WINTER
    if(ii != 2004) track$season[which(track$date >= paste(ii-1,'-12-20',sep='') & track$date < paste(ii-1,'-12-31',sep=''))] <- 4
    track$season[which(track$date >= paste(ii,'-01-01',sep='') & track$date < paste(ii,'-03-20',sep=''))] <- 4
    
    # FALL
    track$season[which(track$date >= paste(ii,'-09-20',sep='') & track$date < paste(ii,'-12-20',sep=''))] <- 3
    
    #SPR
    track$season[which(track$date >= paste(ii,'-03-20',sep='') & track$date < paste(ii,'-06-20',sep=''))] <- 1
    
    #SUM
    track$season[which(track$date >= paste(ii,'-06-20',sep='') & track$date < paste(ii,'-09-20',sep=''))] <- 2
    
  }
  
  spring.idx <- which(track$season == 1)
  summer.idx <- which(track$season == 2)
  fall.idx <- which(track$season == 3)
  winter.idx <- which(track$season == 4)
  
  if(length(spring.idx) > 0) sprRD <- apply(all[,,spring.idx], 1:2, FUN=function(x) sum(x, na.rm=T)) else{sprRD <- apply(all[,,1], 1:2, FUN=function(x) sum(x, na.rm=T)) * NA}
  if(length(summer.idx) > 0) sumRD <- apply(all[,,summer.idx], 1:2, FUN=function(x) sum(x, na.rm=T)) else{sumRD <- apply(all[,,1], 1:2, FUN=function(x) sum(x, na.rm=T)) * NA}
  if(length(fall.idx) > 0) fallRD <- apply(all[,,fall.idx], 1:2, FUN=function(x) sum(x, na.rm=T)) else{fallRD <- apply(all[,,1], 1:2, FUN=function(x) sum(x, na.rm=T)) * NA}
  if(length(winter.idx) > 0) winRD <- apply(all[,,winter.idx], 1:2, FUN=function(x) sum(x, na.rm=T)) else{winRD <- apply(all[,,1], 1:2, FUN=function(x) sum(x, na.rm=T)) * NA}
  
  # sum across time (season) and normalize final surface
  seasonRD <- abind::abind(sprRD, sumRD, fallRD, winRD, along=3)
  seasonRD <- raster::flip(raster::brick(seasonRD, xmn=min(g$lon), xmx=max(g$lon),
                                     ymn=min(g$lat), ymx=max(g$lat), crs=crs, transpose=T), 2)
  
  return(seasonRD=seasonRD)
}