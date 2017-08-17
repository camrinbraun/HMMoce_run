# Compares tracks output from formatTracks.r

compareTracks <- function(df){
  
  # RMSE
  knownVec <- c('spot.lon','spot.lat')
  idx <- which(names(df) %in% knownVec)
  lat.idx <- grep('lat', names(df))
  lat.idx <- lat.idx[-grep('var', names(df)[lat.idx])]
  lat.idx <- lat.idx[!lat.idx %in% idx]
  
  lon.idx <- grep('lon', names(df))
  lon.idx <- lon.idx[-grep('var', names(df)[lon.idx])]
  lon.idx <- lon.idx[!lon.idx %in% idx]
  
  for (i in lat.idx){
    diff.lat <- df[,idx[2]] - df[,i]
    #diff.lat <- abs(diff.lat[-which(is.na(diff.lat))])
    rmse.lat.i <- sqrt(sum((diff.lat ^ 2) / length(df[,idx[2]]), na.rm=T))
    if (i == lat.idx[1]){
      rmse.lat <- rmse.lat.i
    } else{
      rmse.lat <- c(rmse.lat, rmse.lat.i)
    }
  }
  
  for (i in lon.idx){
    diff.lon <- df[,idx[1]] - df[,i]
    rmse.lon.i <- sqrt(sum((diff.lon ^ 2) / length(df[,idx[1]]), na.rm=T))
    if (i == lon.idx[1]){
      rmse.lon <- rmse.lon.i
    } else{
      rmse.lon <- c(rmse.lon, rmse.lon.i)
    }
  }
  
  # MGCE - mean great circle distance using Haversine formula
  # for more info see: https://www.r-bloggers.com/great-circle-distance-calculations-in-r/
  types <- c('ti','tib','kf','kfb','gpe','hmm')
  gcd <- as.data.frame(matrix(NA, ncol=length(types), nrow=nrow(df)))
  names(gcd) <- types
  for (i in 1:length(types)){
    gcd[,i] <- fields::rdist.earth.vec(df[,idx], df[,c(lon.idx[i],lat.idx[i])]) * 1.60934 # results in miles, convert to kms
  }
  
  return(list(rmse.lon = rmse.lon, rmse.lat = rmse.lat, gcd = gcd))
}
