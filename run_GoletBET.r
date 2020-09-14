devtools::load_all('~/work/RCode/HMMoce')
library(raster)
#library(foieGras)
library(tidyverse)

#meta <- read.table('../nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, stringsAsFactors = F)
xlsFiles <- list.files('~/work/Data/Bigeye/WaltTagsAug2020/', full.names = T, recursive = T)
xlsFiles <- xlsFiles[grep('.xls', xlsFiles)]
xlsFiles <- xlsFiles[1]

for (aa in 1:length(xlsFiles)){
  
  xlsfile = xlsFiles[aa]
  locs <- gdata::read.xls(xlsfile, sheet='Lat&Long', skip=1, header=T, stringsAsFactors=F)#[,1:6] # 6 cols
  # dat = MWTextract(tagID = 57508, xlsfile, delta = T, minmax = F)
  day0 = as.POSIXct(locs[3,8], format='%b %d, %Y', tz='UTC') ## tag date
  x0 = as.numeric(locs[5,8:9])#(read_excel(xlsfile, skip = 5, sheet = 'Lat&Long', n_max = 1))[,8:9])
  x0[2] <- x0[2] * -1
  dayT = as.POSIXct(locs[9,8], format='%b %d, %Y', tz='UTC') ## end date
  xT = as.numeric(locs[11,8:9])
  xT[2] <- xT[2] * -1
  
  iniloc <- data.frame(matrix(c(lubridate::day(day0), lubridate::month(day0), lubridate::year(day0), x0[1], x0[2],
                                lubridate::day(dayT), lubridate::month(dayT), lubridate::year(dayT), xT[1], xT[2]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
  tag <- iniloc$date[1]
  pop <- iniloc$date[2]
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  
  ## read MT tag data
  light <- gdata::read.xls(xlsfile, sheet='Sunrise and Sunset Times', skip=1, header=T)[,1:5]
  #lightloc <- gdata::read.xls(xlsfile, sheet='Lat&Long', skip=1, header=T)#[,1:6] # 6 cols
  depth <- gdata::read.xls(xlsfile, sheet='Press Data', skip=1, header=T)[,1:5] # 5 cols
  mmd <- gdata::read.xls(xlsfile, sheet='Press Data (MinMax)', skip=1, header=T)[,1:5] # 5 cols
  temp <- gdata::read.xls(xlsfile, sheet='Temp Data', skip=1, header=T)[,1:5] # 5 cols
  temp_max <- gdata::read.xls(xlsfile, sheet='Temp Data (MinMax)', skip=1, header=T)[,1:5] # 5 cols
  
  ## MMD
  ## create max depth df for calc.bathy
  mmd$Date <- as.POSIXct(mmd$Date.Time, format='%m/%d/%y', tz='UTC')
  mmd$MaxDepth <- mmd$Max.Depth.m. * -1
  mmd <- mmd[,c('Date','MaxDepth')]
  mmd <- mmd[which(mmd$Date <= iniloc$date[2]),]
  
  ## PDT
  ## create depth-temperature profile data for various 3d likelihood calcs
  pdt_dates <- seq.POSIXt(tag, pop, by = 'day') ## generate daily summary of depth-temp profiles
  depth$Date <- as.POSIXct(depth$Date.Time, format='%m/%d/%y %H:%M', tz='UTC')
  depth$Depth <- depth$Depth.m. * -1
  temp$Date <- as.POSIXct(temp$Date.Time, format='%m/%d/%y %H:%M', tz='UTC')
  
  pdt <- merge(depth, temp, by='Date')
  pdt <- pdt[,c('Date','Depth','Temp.C.')]
  names(pdt)[3] <- 'Temperature'
  
  ## generate depth-temp summary from time series
  pdt <- bin_TempTS(pdt, out_dates = pdt_dates, bin_res = 25)
  
  ## SST
  ## assume daily max temp is equivalent to SST
  sst <- temp_max
  sst$Date <- as.POSIXct(sst$Date.Time, format='%m/%d/%y', tz='UTC')
  sst$Temperature <- sst$Max.Temp.C.
  sst <- sst[,c('Date','Temperature')]
  
  ## LIGHT LOCS
  #MWTxy = .getMWTxy(xlsfile, x0, xT, day0, dayT)
  MWTxy = locs[!is.na(locs[, 1]), 1:3]
  MWTxy$Date <- as.POSIXct(MWTxy$Date, format='%b %d, %Y', tz='UTC')
  #dateVec = seq(day0, dayT, by = 'day')
  didx = match(MWTxy$Date, dateVec)
  len = length(dateVec)
  MWTdata = as.data.frame(array(NA, c(len, 3)))
  didx = didx[!is.na(didx)]
  MWTdata[, 1] = dateVec
  MWTdata[, 2:3] = MWTxy[didx, 2:3]
  MWTdata[, 3] = -1 * (MWTdata[, 3])
  MWTdata[1, 2:3] = x0
  MWTdata[len, 2:3] = xT
  names(MWTdata) = c('Date','Latitude','Longitude')
  locs <- MWTdata
  locs$Error.Semi.minor.axis = .7 * 1000 * 111
  #locs$Error.Semi.minor.axis[didx] = fit$var.most.prob.track[,1]*1000*111
  #L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
  
  ## IT IS HIGHLY RECOMMENDED YOU CHECK AND MANUALLY FILTER OUTPUT LOCATIONS FROM MICROWAVE TAGS!!
  ## THESE ARE VERY RAW LIGHT-BASED ESTIMATES AND WHEN OUTLIERS ARE INCLUDED THEY CAN DRAMATICALLY IMPACT HMMOCE OUTPUTS
  ## here we started with 37 estimates and filtered to 8 reasonably useful ones
  plot(locs$Date, locs$Longitude, pch=16); world(add=T)
  if (aa == 1) locs <- locs[which(locs$Longitude < -40 & locs$Longitude > -80),]
  plot(locs$Longitude, locs$Latitude); world(add=T)
  
  #invisible(readline(prompt=paste('Check plots. If passes QC, press [enter] to continue.', sep='')))
  
  
  if (aa == 1){
    sp.lim <- list(lonmin = -80,
                   lonmax = -60,
                   latmin = 15,
                   latmax = 45)
    
  } else if (aa == 2){
    sp.lim <- list(lonmin = -80,
                   lonmax = -60,
                   latmin = 15,
                   latmax = 45)
    
  } else if (aa == 3){
    sp.lim <- list(lonmin = -80,
                   lonmax = -60,
                   latmin = 15,
                   latmax = 45)
  }
    
  ## setup the spatial grid to base likelihoods on
  locs.grid <- setup.locs.grid(sp.lim, res='quarter')
  
  udates <- seq.Date(as.Date(tag), as.Date(pop), by = 'day')
  
  setwd('~/work/RCode/HMMoce_run/57508/')
  dir <- getwd()
  #load('./141259_tryHMMoce_20200724.rda')
  #load('./141259_L.res_20200727.rda')
  
  sst.dir <- paste0(dir, '/EnvData/sst/')
  if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
  for (i in 1:length(udates)) get.env(udates[i], filename='oisst', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
  
  hycom.dir <- paste0(dir,'/EnvData/hycom/')
  if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
  for (i in 1:length(udates)) get.env(udates[i], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
  
  bathy.dir <- paste0(dir, '/EnvData/bathy/')
  if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
  #if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
  bathy <- HMMoce::get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir, res=1)
  #} else{ ## OR (once downloaded and reading the .nc later)
  bathy <- irregular_ncToRaster(paste0(bathy.dir, 'bathy.nc'), varid = 'topo')
  #}
  
  
  L.light <- calc.lightloc(locs, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE)
  
  L.sst <- calc.sst(sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1, focalDim = 3)
  
  
  # OCEAN HEAT CONTENT (INTEGRATED PDTS)
  L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  #L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  
  # HYCOM PROFILE BASED LIKELIHOODS
  L.hycom.se <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
  #L.hycom.par <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:5], use.se = F, ncores=2)
  
  bathy_resamp <- raster::resample(bathy, L.sst) # or whatever grid makes sense to resample to
  L.bathy <- calc.bathy(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')
  
  ## make list of rasters
  L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.hycom.se = L.hycom.se, L.bathy = L.bathy)
  
  
}
