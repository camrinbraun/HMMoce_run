dataDir <- '~/ebs/Data/BaskingSharks/batch/'
library(raster); library(rgeos); library(RNetCDF)
meta <- read.table(paste(dataDir, 'bask_metadata_v2.csv',sep=''), sep=',', header=T)
meta$PopDate <- as.Date(meta$PopDate, format='%m/%d/%y')

profs <- list()
for (tt in 35:nrow(meta)){
  ptt <- meta$PTT[tt]
  dte <- meta$PopDate[tt]
  if (meta$bnds[tt] == 'big'){
    ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/big/bask_big_', dte, '.nc', sep='')
  } else{
    ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/small/bask_', dte, '.nc', sep='')
  }
  
  
  if(!file.exists(ncName)){
    if (meta$bnds[tt] == 'big'){
      ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/big/bask_big_', dte-1, '.nc', sep='')
    } else{
      ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/small/bask_', dte-1, '.nc', sep='')
    }
  }
  
  # open nc and get the indices for the vars
  nc1 =  RNetCDF::open.nc(ncName)
  ncnames = NULL
  nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
  for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
  temp.idx <- grep('temp', ncnames, ignore.case=TRUE) - 1
  lat.idx <- grep('lat', ncnames, ignore.case=TRUE) - 1
  lon.idx <- grep('lon', ncnames, ignore.case=TRUE) - 1
  dep.idx <- grep('dep', ncnames, ignore.case=TRUE) - 1
  
  # get attributes, if they exist
  ncatts <- NULL
  nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
  for(ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1, temp.idx, ii)$name
  scale.idx <- grep('scale', ncatts, ignore.case=TRUE) - 1
  if(length(scale.idx) != 0){
    scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=scale.idx)
  } else{
    scale <- 1
  }
  off.idx <- grep('off', ncatts, ignore.case=TRUE) - 1
  if(length(off.idx) != 0){
    offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=off.idx)
  } else{
    offset <- 1
  }
  
  # get and check the vars
  depth <- RNetCDF::var.get.nc(nc1, dep.idx)
  lon <- RNetCDF::var.get.nc(nc1, lon.idx)
  if(length(dim(lon)) == 2) lon <- lon[,1]
  if(!any(lon < 180)) lon <- lon - 360
  lat <- RNetCDF::var.get.nc(nc1, lat.idx)
  if(length(dim(lat)) == 2) lat <- lat[1,]
  
  # position indices
  elo <- which.min(abs(lon - meta$PopLong[tt]))
  ela <- which.min(abs(lat - meta$PopLat[tt]))
  
  # get profile
  temp <- var.get.nc(nc1, 'water_temp', start=c(elo,ela,1,1), count=c(1,1,NA,NA)) * scale + offset
  
  profs[[tt]] <- list(ptt=ptt, temp=temp, depth=depth, pos=c(meta$PopLong, meta$PopLat))
}


for (ii in 1:nrow(meta)){
  
  dataDir <- '~/ebs/Data/BaskingSharks/batch/'
  ptt <- meta$PTT[ii]
  setwd(paste(dataDir, '/', ptt, sep=''))
  
  iniloc <- data.frame(matrix(c(meta$TagDay[ii], meta$TagMonth[ii], meta$TagYear[ii], meta$TagLat[ii], meta$TagLong[ii], 
                                meta$PopDay[ii], meta$PopMonth[ii], meta$PopYear[ii], meta$PopLat[ii], meta$PopLong[ii]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  
  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day')) 
  
  # READ IN DATA FROM WC FILES
  myDir <- paste(dataDir, ptt, '/', sep='')
  setwd(myDir)

  # depth-temp profile data
  pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
  pdt.udates <- pdt$udates; pdt <- pdt$data
  tail(pdt)
  print(meta$PopDate[ii])
  pdt$Date <- as.Date(pdt$Date, format = findDateFormat(pdt$Date))
  
  pdt.i <- pdt[which(pdt$Date %in% as.Date('2012-05-30')),c(2,5:7)]
  
  profs[[ii]]$pdt.i <- pdt.i

#}
ii=ii+1


#> setwd('~/ebs/Data/BaskingSharks/batch/')
#> save(profs, file='bask_profile_list.rda')
  