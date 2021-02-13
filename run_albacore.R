#================
## INTRO STUFF
#================

#devtools::load_all('../HMMoce')
devtools::install_github('camrinbraun/HMMoce', ref='dev')
library(HMMoce)
#library(raster)
library(tidyverse)
#devtools::install_github('camrinbraun/tags2etuff')
#library(tags2etuff)
devtools::load_all('../tags2etuff') ## for etuff functionality
#devtools::load_all('../analyzePSAT')
source('../analyzePSAT/R/make360.R')
library(fields) ## for quick mapping
library(raster)
<<<<<<< HEAD
library(foreach)
#library(GA) ## for genetic algorithm

# aws s3 sync /home/rstudio/ebs/RCode/HMMoce_run/ s3://braunmpg/RCode/HMMoce_run/ --exclude "*" --include "*_res*"
# aws s3 sync /home/rstudio/ebs/RCode/HMMoce_run/ s3://braunmpg/RCode/HMMoce_run/ --exclude "*" --include "*L_tt*"
# aws s3 sync /home/rstudio/ebs/RCode/HMMoce_run/ s3://braunmpg/RCode/HMMoce_run/ --exclude "*" --include "*L.res*"
# aws s3 sync /home/rstudio/ebs/RCode/HMMoce_run/ s3://braunmpg/RCode/HMMoce_run/ --exclude "*" --include "*L.rasters*"
=======
#library(GA) ## for genetic algorithm
>>>>>>> 2597194652866eb91bd5bf3a16df89994df86710

base_dir <- '~/ebs/RCode/HMMoce_run'
setwd(base_dir)

meta <- read.table('../nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, stringsAsFactors = F)
meta <- meta %>% filter(platform == 'Thunnus alalunga' & etuff == 1)
meta$time_coverage_start <- as.POSIXct(meta$time_coverage_start, tz = 'UTC')
meta$time_coverage_end <- as.POSIXct(meta$time_coverage_end, tz = 'UTC')

## a set of manually defined spatial bounds for each fish
limits <- read.table('~/ebs/Data/albacore/alb_limits.csv', sep=',', header=T)

## which light-based lon estimates need to be removed from further analysis?
drop_light <- readRDS('~/ebs/Data/albacore/drop_light_CDB.rds')
<<<<<<< HEAD


## redo SST, calc new L.sst done
id_idx <- which(meta$instrument_name %in%
                  c("172419_2011_1090251",
                    "172419_2003_390191"))



for (i in 4:nrow(meta)){

#for (i in 10:nrow(meta)){
=======

#for (i in 1:nrow(meta)){

for (i in 3:nrow(meta)){
>>>>>>> 2597194652866eb91bd5bf3a16df89994df86710
  
  if (class(drop_light[[i]]) == 'logical'){
    drop_vec <- NA
  } else if (class(drop_light[[i]]) == 'list'){
    drop_vec <- c(drop_light[[i]]$by_hand)
  } else if (class(drop_light[[i]]) == 'integer'){
    drop_vec <- drop_light[[i]]
  } else{
    stop('drop_light failed')
  }
  if (length(drop_vec) == 0) drop_vec <- NA
  
  ## temporal bounds
  iniloc <- data.frame(matrix(c(lubridate::day(meta$time_coverage_start[i]),
                                lubridate::month(meta$time_coverage_start[i]),
                                lubridate::year(meta$time_coverage_start[i]),
                                meta$geospatial_lat_start[i],
                                make360(meta$geospatial_lon_start[i]),
                                lubridate::day(meta$time_coverage_end[i]),
                                lubridate::month(meta$time_coverage_end[i]),
                                lubridate::year(meta$time_coverage_end[i]),
                                meta$geospatial_lat_end[i],
                                make360(meta$geospatial_lon_end[i])),
                              nrow = 2, ncol = 5, byrow = T))
  
  names(iniloc) <- list('day','month','year','lat','lon')
  iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
  tag <- iniloc$date[1]
  pop <- iniloc$date[2]
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  
  sp.lim <- list(lonmin = make360(limits$lonmin[which(limits$instrument_name == meta$instrument_name[i])]),
                 lonmax = make360(limits$lonmax[which(limits$instrument_name == meta$instrument_name[i])]),
                 latmin = limits$latmin[which(limits$instrument_name == meta$instrument_name[i])],
                 latmax = limits$latmax[which(limits$instrument_name == meta$instrument_name[i])])
  
  ## setup the spatial grid to base likelihoods on
  locs.grid <- setup.locs.grid(sp.lim, res='quarter')
  
  work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
  if (!dir.exists(work_dir)) dir.create(work_dir)
  setwd(work_dir)
  
  sst.dir <- paste0('./tmp/sst/')
  #if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
  #for (b in 1:length(dateVec)){
  #  if (!file.exists(paste0(sst.dir, '/', 'oi_', dateVec[b], '.nc'))) get.env(dateVec[b], filename='oi', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
  #}
  if (length(list.files(sst.dir)) != length(dateVec)){
    system(paste0('aws s3 cp s3://braunmpg/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/sst/ /home/rstudio/ebs/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/sst/'))
    if (length(list.files(sst.dir)) != length(dateVec)) stop('Check that all daily sst data is available for this individual run.')
  }
  
  hycom.dir <- paste0('./tmp/hycom/')
  if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
<<<<<<< HEAD
  ## check for hycom. if its not available pull it down from the s3 bucket
  if (length(list.files(hycom.dir)) != length(dateVec)){
    #system(
      paste0('aws s3 sync s3://braunmpg/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/ /home/rstudio/ebs/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/')
    #)
    if (length(list.files(hycom.dir)) != length(dateVec)) stop('Check that all daily hycom data is available for this individual run.')
  }
  
  ## this is how we get hycom. it will skip days for which data already exists
  for (b in 1:length(dateVec)){
    setwd(work_dir)
=======
  
  ## check for hycom. if its not available pull it down from the s3 bucket
  #if (length(list.files(hycom.dir)) != length(dateVec)){
  #  system(paste0('aws s3 cp s3://braun-data/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/ /home/rstudio/ebs/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/'))
  #  if (length(list.files(hycom.dir)) != length(dateVec)) stop('Check that all daily hycom data is available for this individual run.')
  #}
  
  ## this is how we get hycom. it will skip days for which data already exists
  for (b in 1:length(dateVec)){
>>>>>>> 2597194652866eb91bd5bf3a16df89994df86710
    if (!file.exists(paste0(hycom.dir, '/', 'hycom_', dateVec[b], '.nc'))) {
      print(b)
      try(
<<<<<<< HEAD
      get.env(dateVec[b], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir), TRUE)
      
  #    #if (b %% 10 == 0){
  #      #removeTmpFiles(); 
  #    removeTmpFiles(h=0)
  #    file.remove(list.files(tempdir(), recursive=TRUE))
  #    #gc()
      }
    }
  #}
    
  #bathy.dir <- paste0('./tmp/bathy/')
  #if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
  #if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
  #  bathy <- HMMoce::get.bath.data(sp.lim, save.dir = bathy.dir, res=1)
  #} else{ ## OR (once downloaded and reading the .nc later)
    bathy <- raster::raster('~/ebs/EnvData/bathy/global_bathy_0.01.nc')
    bathy <- raster::crop(bathy, raster::extent(unlist(sp.lim)))
  #  raster::writeRaster(bathy, paste0(bathy.dir, '/bathy.nc'), format='CDF', overwrite=TRUE, varname = 'topo')                
  #}
  #print(i); gc()

  data_dir <- paste0('~/ebs/Data/albacore/', meta$instrument_name[i], '/cdb/')
  #  data_dir <- paste0('~/Google Drive File Stream/My Drive/Albacore - All Data/data/', meta$instrument_name[i], '/cdb/')
  etuff_file <- paste(data_dir, meta$instrument_name[i], '_eTUFF.txt', sep='')
  
  etuff <- read_archival(etuff_file)
  
  ## sst is Date, Temperature
  ## most tags have some reported SST metric but easier to create our own with archival tags
  #sstVars <- varNames[grep('sst', varNames)]
  #sst1 <- archival_to_etuff(etuff$etuff, vars = c('DateTime', sstVars))
  series <- get_series(etuff)
  #rm(etuff)
  sst <- data.frame(series %>% filter(!is.na(temperature)) %>% 
                      group_by(as.Date(DateTime)) %>%
                      summarise(n=n(),
                                length_temp = length(temperature[which.min(depth)]),
                                which_min = which.min(depth),
                                sst = temperature[which.min(depth)],
                                sst_depth = depth[which.min(depth)],
                                sst_maxT = max(temperature, na.rm=T),
                                sst_mean = mean(temperature[which(depth <= 2)])))
  sst <- sst %>% dplyr::select(as.Date.DateTime., sst, sst_depth)
  names(sst) <- c('Date','Temperature', 'sst_depth')
  sst$Date <- as.POSIXct(sst$Date, tz='UTC')
   
  
  sst2 <- archival_to_etuff(etuff$etuff, vars ='sst')
  sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T), mean_sst = mean(sstMean, na.rm=T))
  #sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T))
  sst2$Date <- as.POSIXct(sst2$`as.Date(DateTime)`, tz='UTC')
  sst2 <- data.frame(sst2)
  for (b in 1:nrow(sst2)){
    if(is.nan(sst2$Temperature[b]) & !is.nan(sst2$mean_sst[b])) sst2$Temperature[b] <- sst2$mean_sst[b]
  }
  sst2 <- sst2 %>% dplyr::select(Date, Temperature)
  
  comp_sst <- merge(sst, sst2, by='Date')
  ggplot(comp_sst, aes(x=Temperature.x, y=Temperature.y, colour=sst_depth)) + geom_point() +
    geom_abline(slope=1, intercept=0)
  
  plot(comp_sst$Temperature.x, comp_sst$Temperature.y, xlab='eTuff SST', ylab='lotek sst')
  abline(a=0, b=1)
  
  
  wc_sst <- data.table::fread('~/ebs/Data/albacore/172419_2011_1090269/cdb/1090269-SST.csv', sep=',')
  wc_sst$Date <- as.POSIXct(wc_sst$Date, format='%H:%M:%S %d-%b-%Y', tz='UTC')
  wc_sst <- wc_sst %>% group_by(as.Date(Date)) %>% summarise(n=n(), mean_sst = mean(Temperature), min_sst = min(Temperature), max_sst = max(Temperature))
  sst <- sst %>% group_by(as.Date(Date)) %>% summarise(n=n(), mean_sst = mean(Temperature), min_sst = min(Temperature), max_sst = max(Temperature))
  sst2 <- archival_to_etuff(etuff$etuff, vars =)
  sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), mean_sst = mean(sst, na.rm=T), min_sst = min(sst, na.rm=T), max_sst = max(sst, na.rm=T))
  names(sst)[1] <- 'Date'; names(sst2)[1] <- 'Date'; names(wc_sst)[1] <- 'Date'
  
  wc_series <- merge(wc_sst, sst, by='Date')
  wc_etuff <- merge(wc_sst, sst2, by='Date')
  
  plot(wc_series$mean_sst.x, wc_series$mean_sst.y, xlab='WC SST', ylab='eTUFF series')
  abline(a=0, b=1)
  
  plot(wc_etuff$mean_sst.x, wc_etuff$mean_sst.y, xlab='WC SST', ylab='eTUFF sst')
  abline(a=0, b=1)
  
  ## mmd is Date, MaxDepth
  mmd <- data.frame(series %>% filter(!is.na(temperature)) %>% 
                      group_by(as.Date(DateTime)) %>%
                      summarise(MaxDepth = max(depth)))
  names(mmd) <- c('Date','MaxDepth')
  mmd$MaxDepth <- round(mmd$MaxDepth, 0)
  mmd$Date <- as.POSIXct(mmd$Date, tz='UTC')
  
  ## depth-temp profiles: Date, Depth, MeanTemp
  pdt <- data.frame(series %>% filter(!is.na(temperature)))
  pdt$Date <- as.Date(pdt$DateTime)
  hycom_depth <- c(0, 2, 4, 6, 8, 10, 12, 15, 20, 25,
                   30, 35, 40, 45, 50, 60, 70, 80, 90,
                   100, 125, 150, 200, 250, 300, 350, 
                   400, 500, 600, 700, 800, 900, 1000,
                   1250, 1500, 2000, 2500, 3000, 4000, 5000)
  pdt$depth[which(pdt$depth < 0)] <- 0
  pdt$depthInterval <- findInterval(pdt$depth, hycom_depth)
  pdt <- data.frame(pdt %>% group_by(Date, depthInterval) %>% 
                      summarise(MinTemp = min(temperature), MaxTemp = max(temperature), 
                                MeanTemp = mean(temperature), n=n()))
  pdt$depth <- hycom_depth[pdt$depthInterval]
  pdt <- pdt[,c('Date','depth','MinTemp','MaxTemp')]
  names(pdt)[2] <- 'Depth'
  pdt$MinTemp <- round(pdt$MinTemp, 1)
  pdt$MaxTemp <- round(pdt$MaxTemp, 1)
  pdt$Date <- as.POSIXct(pdt$Date, tz='UTC')
  
  ## get light locs
  tr <- get_track(etuff)
  tr$lon2 <- make360(tr$longitude)
  auto_drop <- which(tr$lon2 < 125 | tr$lon2 > 260)
  if (length(auto_drop) > 0) tr <- tr[-auto_drop,]
  if (all(!is.na(drop_vec))) tr <- tr[-drop_vec,]
  if (meta$manufacturer[i] == 'Lotek'){
    ## then we won't have error estimates. Musyl error estimates will be applied automatically in the likelihood calc
    tr <- tr[,c('DateTime','lon2','latitude')]
    names(tr) <- c('Date','Longitude','Latitude')
  } else{
    ## its WC and we do have them
    tr <- tr[,c('DateTime','lon2','latitude','longitudeError')]
    names(tr) <- c('Date','Longitude','Latitude','longitudeError')
  }
  if (i == 8) tr <- tr %>% filter(Date != as.POSIXct('2013-06-25', tz='UTC')) ## weird one we missed in filtering
  
  world <- map_data('world2')
  xl <- c(min(tr$Longitude) - 5, max(tr$Longitude) + 5)
  yl <- c(min(tr$Latitude) - 5, max(tr$Latitude) + 5)
  
  ## simple map of move data just to double-check before calculating light-based likelihoods
  ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
    coord_fixed(xlim=xl, ylim=yl) + xlab('') + ylab('') +
    geom_point(data=tr, aes(x=Longitude, y=Latitude, colour=Date)) +
    geom_point(data = meta[i,], aes(x=make360(geospatial_lon_start), y=geospatial_lat_start), col='green') +
    geom_point(data = meta[i,], aes(x=make360(geospatial_lon_end), y=geospatial_lat_end), col='red') #+
    #geom_point(data = res$tr, aes(x=lon, y=lat))
    
  
  ## LIGHT
  L.light <- calc.lightloc(tr, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE)
  
  ## SST
  #L.sst <- calc.sst.par(sst, filename='oi', sst.dir = sst.dir, dateVec = dateVec[1:10], sens.err = 3, focalDim = 3)
  L.sst <- calc.sst.par(sst2, filename='oi', sst.dir = sst.dir, dateVec = dateVec, sens.err = 3, focalDim = 3)
  ## ~9 mins for Lotek 0394 using 8 cores (m4.2xlarge)
  ## 23 mins for Lotek 0396 using 8 cores (m4.2xlarge)
  ## 5 mins for Lotek 0396 using 40 cores (m4.10xlarge)
  
  load('172419_2003_390191_L.res_20210114.rda')
  
  pdf('compare_sst.pdf', height=10, width=12)
  par(mfrow=c(2,1))
  for (b in 1:length(dateVec)){
    plot(L.res$L.rasters$L.sst[[b]], main='original sst')
    world(add=T, wrap=c(0,360))
    plot(L.sst[[b]], main='new sst')
    world(add=T, wrap=c(0,360))
  }
 
  dev.off()
  
  # OCEAN HEAT CONTENT (INTEGRATED PDTS)
  L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec[3], isotherm = '', use.se = F)
  #L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, 
                        isotherm = '', use.se = F)#, ncores = 30)
  L.ohc
  L.ohc <- raster::shift(L.ohc, dx=360)
  #writeRaster(L.ohc, filename = paste0(work_dir, 'L.ohc.grd'), format='raster')
  ## 50 mins for Lotek 0396 using 40 cores (m4.10xlarge)

  # WORLD OCEAN ATLAS-BASED LIKELIHOODS
  #L.woa <- calc.woa(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
  #L.woa <- calc.woa.par(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
  
  # HYCOM PROFILE BASED LIKELIHOODS
  #L.hycom <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:3], use.se = T)
  #L.hycom.se <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:3], use.se = T)
  #L.hycom.par <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:5], use.se = F, ncores=2)
=======
        get.env(dateVec[b], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir), TRUE)
      #if (b %% 10 == 0){
      #removeTmpFiles(); 
      removeTmpFiles(h=0)
      file.remove(list.files(tempdir(), recursive=TRUE))
      #gc()
      #}
    }
  }
}

bathy.dir <- paste0('./tmp/bathy/')
if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
#if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
#  bathy <- HMMoce::get.bath.data(sp.lim, save.dir = bathy.dir, res=1)
#} else{ ## OR (once downloaded and reading the .nc later)
bathy <- raster::raster('~/ebs/EnvData/bathy/global_bathy_0.01.nc')
bathy <- raster::crop(bathy, raster::extent(unlist(sp.lim)))
#  raster::writeRaster(bathy, paste0(bathy.dir, '/bathy.nc'), format='CDF', overwrite=TRUE, varname = 'topo')                
#}
#print(i); gc()




data_dir <- paste0('~/ebs/Data/albacore/', meta$instrument_name[i], '/cdb/')
#  data_dir <- paste0('~/Google Drive File Stream/My Drive/Albacore - All Data/data/', meta$instrument_name[i], '/cdb/')
etuff_file <- paste(data_dir, meta$instrument_name[i], '_eTUFF.txt', sep='')

etuff <- read_archival(etuff_file)

## sst is Date, Temperature
## most tags have some reported SST metric but easier to create our own with archival tags
#sstVars <- varNames[grep('sst', varNames)]
#sst1 <- archival_to_etuff(etuff$etuff, vars = c('DateTime', sstVars))
series <- get_series(etuff)
#rm(etuff)
sst <- data.frame(series %>% filter(!is.na(temperature)) %>% 
                    group_by(as.Date(DateTime)) %>%
                    summarise(sst = temperature[which.min(depth)]))
names(sst) <- c('Date','Temperature')
sst$Date <- as.POSIXct(sst$Date, tz='UTC')

## mmd is Date, MaxDepth
mmd <- data.frame(series %>% filter(!is.na(temperature)) %>% 
                    group_by(as.Date(DateTime)) %>%
                    summarise(MaxDepth = max(depth)))
names(mmd) <- c('Date','MaxDepth')
mmd$MaxDepth <- round(mmd$MaxDepth, 0)
mmd$Date <- as.POSIXct(mmd$Date, tz='UTC')

## depth-temp profiles: Date, Depth, MeanTemp
pdt <- data.frame(series %>% filter(!is.na(temperature)))
pdt$Date <- as.Date(pdt$DateTime)
hycom_depth <- c(0, 2, 4, 6, 8, 10, 12, 15, 20, 25,
                 30, 35, 40, 45, 50, 60, 70, 80, 90,
                 100, 125, 150, 200, 250, 300, 350, 
                 400, 500, 600, 700, 800, 900, 1000,
                 1250, 1500, 2000, 2500, 3000, 4000, 5000)
pdt$depth[which(pdt$depth < 0)] <- 0
pdt$depthInterval <- findInterval(pdt$depth, hycom_depth)
pdt <- data.frame(pdt %>% group_by(Date, depthInterval) %>% 
                    summarise(MinTemp = min(temperature), MaxTemp = max(temperature), 
                              MeanTemp = mean(temperature), n=n()))
pdt$depth <- hycom_depth[pdt$depthInterval]
pdt <- pdt[,c('Date','depth','MinTemp','MaxTemp')]
names(pdt)[2] <- 'Depth'
pdt$MinTemp <- round(pdt$MinTemp, 1)
pdt$MaxTemp <- round(pdt$MaxTemp, 1)
pdt$Date <- as.POSIXct(pdt$Date, tz='UTC')

## get light locs
tr <- get_track(etuff)
tr$lon2 <- make360(tr$longitude)
auto_drop <- which(tr$lon2 < 125 | tr$lon2 > 260)
if (length(auto_drop) > 0) tr <- tr[-auto_drop,]
if (all(!is.na(drop_vec))) tr <- tr[-drop_vec,]
if (meta$manufacturer[i] == 'Lotek'){
  ## then we won't have error estimates. Musyl error estimates will be applied automatically in the likelihood calc
  tr <- tr[,c('DateTime','lon2','latitude')]
  names(tr) <- c('Date','Longitude','Latitude')
} else{
  ## its WC and we do have them
  tr <- tr[,c('DateTime','lon2','latitude','longitudeError')]
  names(tr) <- c('Date','Longitude','Latitude','longitudeError')
}

world <- map_data('world2')
xl <- c(min(tr$Longitude) - 5, max(tr$Longitude) + 5)
yl <- c(min(tr$Latitude) - 5, max(tr$Latitude) + 5)

## simple map of move data just to double-check before calculating light-based likelihoods
ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
  coord_fixed(xlim=xl, ylim=yl) + xlab('') + ylab('') +
  geom_point(data=tr, aes(x=Longitude, y=Latitude, colour=Date)) +
  geom_point(data = meta[i,], aes(x=make360(geospatial_lon_start), y=geospatial_lat_start), col='green') +
  geom_point(data = meta[i,], aes(x=make360(geospatial_lon_end), y=geospatial_lat_end), col='red') +
  geom_point(data = res$tr, aes(x=lon, y=lat))


## LIGHT
L.light <- calc.lightloc(tr, locs.grid = locs.grid, dateVec = dateVec, errEll = TRUE)

## SST
L.sst <- calc.sst.par(sst, filename='oi', sst.dir = sst.dir, dateVec = dateVec, sens.err = 3, focalDim = 3)
## ~9 mins for Lotek 0394 using 8 cores (m4.2xlarge)

# OCEAN HEAT CONTENT (INTEGRATED PDTS)
#L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
## task 1 failed - "invalid type (NULL) for variable 'pdt.i$mintemp'"

# WORLD OCEAN ATLAS-BASED LIKELIHOODS
#L.woa <- calc.woa(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
#L.woa <- calc.woa.par(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)

# HYCOM PROFILE BASED LIKELIHOODS
#L.hycom <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:3], use.se = T)
#L.hycom.se <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:3], use.se = T)
#L.hycom.par <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:5], use.se = F, ncores=2)

#L.bt <- calc.bottomTemp(tag.bt, dateVec[1:5], focalDim = 3, sens.err = 1, bt.dir = sst.dir, filename = 'oisst', varName = 'sst')


## lets see how the below works using the full-res bathymetry calculating likelihoods in parallel
## you can see i setup an if() to coarsen the bathy grid for those that moved offshore but
## i changed my mind given that all fish spend some time coastal and those particular movements
## likely will benefit from higher-res bathymetry
#if (sp.lim$lonmin < 210){
## if the fish made extensive movements offshore (to the west), coarsen the grid to ease computation
#  bathy_resamp <- raster::resample(bathy, L.sst) # or whatever grid makes sense to resample to
#  L.bathy <- calc.bathy(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')
#} else{
## otherwise use the full resolution bathymetry for a fish that stayed more coastal
L.bathy <- calc.bathy.par(mmd, bathy, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')
#}

## make list of rasters
L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.bathy = L.bathy)

## resample rasters
resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
save(L.res, file=paste0(meta$instrument_name[i],'_L.res_', format(Sys.Date(), '%Y%m%d'), '.rda'))
save(L.rasters, file=paste0(meta$instrument_name[i],'_L.rasters_', format(Sys.Date(), '%Y%m%d'), '.rda'))

load(file=paste0(meta$instrument_name[i],'_L.res_', format(Sys.Date(), '%Y%m%d'), '.rda'))


run_idx <- list(c(1:4),
                c(1,2,4),
                c(1,3,4))

## each run idx combination
for (tt in 1:length(run_idx)){
  
  #iniloc$lon <- make360(iniloc$lon)
  L <- make.L(L.res$L.rasters[run_idx[[tt]]], iniloc, dateVec)
  #L.mle <- coarse.L(L, L.res$L.rasters)$L.mle
  #g.mle <- coarse.L(L, L.res$L.rasters)$g.mle
  
  #pars.optim <- opt.params(pars.init = c(2,.2,.6,.8), 
  #                         lower.bounds = c(0.1, 0.001, .1, .1),  
  #                         upper.bounds = c(6, .6, .9, .9), 
  #                         g = L.res$g, 
  #                         #g = g.mle, 
  #                         L = L, 
  #                         #L = L.mle, 
  #                         alg.opt = 'optim', 
  #                         write.results = FALSE)
  
  #pars.optim.mle <- opt.params(pars.init = c(2,.2,.6,.8), 
  #                             lower.bounds = c(0.1, 0.001, .1, .1),  
  #                             upper.bounds = c(6, .6, .9, .9), 
  #                             #g = L.res$g, 
  #                             g = g.mle, 
  #                             #L = L, 
  #                             L = L.mle, 
  #                             alg.opt = 'optim', 
  #                             write.results = FALSE)
  
  ## about 22 mins per iteration on blue shark 141259
  ## about 1.5 mins with MLE grid
  ## final vals way different between the different grids
  
  #pars.nlminb <- opt.params(pars.init = c(2,.2,.6,.8), 
  #                          lower.bounds = c(0.1, 0.001, .1, .1), 
  #                          upper.bounds = c(5, .5, .9, .9), 
  #                          g = L.res$g, 
  #                          #g = g.mle, 
  #                          L = L, 
  #                          #L = L.mle, 
  #                          alg.opt = 'nlminb', 
  #                          write.results = FALSE)
  ## about 30 mins per iteration on blue shark 141259
  
  #ceiling(parallel::detectCores() * .9)
  #pars.ga.one <- opt.params(pars.init = c(2), 
  #                          lower.bounds = c(1), 
  #                          upper.bounds = c(8), 
  #                          g = L.res$g, 
  #                          #g = g.mle, 
  #                          L = L, 
  #                          #L = L.mle, 
  #                          alg.opt = 'ga', 
  #                          write.results = FALSE,
  #                          ncores = ceiling(parallel::detectCores() * .9))
  #ncores = 2)
  
  #pars.ga.mle <- opt.params(pars.init = c(2,.2,.6,.8), 
  #                          lower.bounds = c(0.1, 0.001, .1, .1), 
  #                          upper.bounds = c(6, .6, .9, .9), 
  #                          #g = L.res$g, 
  #                          g = g.mle, 
  #                          #L = L, 
  #                          L = L.mle, 
  #                          alg.opt = 'ga', 
  #                          write.results = FALSE,
  #                          ncores = ceiling(parallel::detectCores() * .9))
  #ncores = 2)
  
  pars.ga <- HMMoce:::opt.params(pars.init = c(2,.2,.6,.8), 
                                 lower.bounds = c(0.1, 0.001, .1, .1), 
                                 upper.bounds = c(6, .6, .9, .9),
                                 max_iter = 10,
                                 run = 10,
                                 g = L.res$g, 
                                 #g = g.mle, 
                                 L = L, 
                                 #L = L.mle, 
                                 alg.opt = 'ga', 
                                 write.results = FALSE,
                                 ncores = 8)
  #ncores = ceiling(parallel::detectCores() * .9))
  pars <- pars.ga$par
  sigmas = pars[1:2]
  sizes = rep(ceiling(sigmas[1]*4),2)
  pb = pars[3:4]
  muadvs = c(0,0)
  
  ## behav 1
  if(sizes[1]%%2==0){sizes[1]=sizes[1]+1}
  ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
  while(ss<.999){
    sizes[1]=sizes[1]+2
    ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
  }
  K1 <- HMMoce:::gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
  rm(ss)
  K1 <- HMMoce:::mask.K(K1)
  
  ## behav 2
  if(sizes[2] %% 2 == 0) sizes[2] <- sizes[2] + 1
  ss <- sum(gausskern.nostd(sizes[2], sigmas[2], muadv = muadvs[2]))
  while(ss < .999){
    sizes[2] <- sizes[2] + 2
    ss <- sum(gausskern.nostd(sizes[2], sigmas[2], muadv = muadvs[2]))
  }
  K2 <- HMMoce:::gausskern.pg(sizes[2], sigmas[2], muadv = muadvs[2])
  K2 <- HMMoce:::mask.K(K2)
>>>>>>> 2597194652866eb91bd5bf3a16df89994df86710
  
  P <- matrix(c(pb[1], 1 - pb[1], 1 - pb[2], pb[2]), 2, 2, byrow = TRUE)
  
<<<<<<< HEAD
  
  ## lets see how the below works using the full-res bathymetry calculating likelihoods in parallel
  ## you can see i setup an if() to coarsen the bathy grid for those that moved offshore but
  ## i changed my mind given that all fish spend some time coastal and those particular movements
  ## likely will benefit from higher-res bathymetry
  #if (sp.lim$lonmin < 210){
    ## if the fish made extensive movements offshore (to the west), coarsen the grid to ease computation
    # bathy_resamp <- raster::resample(bathy, L.sst) # or whatever grid makes sense to resample to
  #  L.bathy <- calc.bathy(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')
  #} else{
    ## otherwise use the full resolution bathymetry for a fish that stayed more coastal
    L.bathy <- calc.bathy.par(mmd, bathy, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')
  #}
  
  gc()
    
  #L.ohc <- raster::brick('./L.ohc.grd')
  ## make list of rasters
  L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.bathy = L.bathy)
  save(L.rasters, file=paste0(meta$instrument_name[i],'_L.rasters_', format(Sys.Date(), '%Y%m%d'), '.rda'))
    
  ## resample rasters
  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  #L.rasters$L.ohc <- raster::shift(L.rasters$L.ohc, dx = 360)
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  L.res$L.rasters$L.sst <- raster::resample(L.sst, L.res$L.rasters$L.sst)
  save(L.res, file=paste0(meta$instrument_name[i],'_L.res_', format(Sys.Date(), '%Y%m%d'), '.rda'))

}






id_idx <- which(meta$instrument_name %in%
                  c("172419_2004_B2381",
                  "172419_2006_D1464",
                  "172419_2003_A1973",
                  "172419_2004_A2088"))

## redo SST, needs new L.sst
id_idx <- which(meta$instrument_name %in%
                  c("172419_2011_1190241",
                    "172419_2015_1490108",
                    "172419_2004_A2088",
                    "172419_2004_B2605",
                    "172419_2003_390167",
                    "172419_2003_390173"))

#for (i in 24:nrow(meta)){
for (i in id_idx){

  #if (meta$instrument_name[i] %in% c("172419_2004_B2381",
  #                                   "172419_2006_D1464",
  #                                   "172419_2003_A1973",
  #                                   "172419_2004_A2088",
  #                                   "172419_2003_390191")) next ## these need cut dates or other special attention
   
  ## temporal bounds
  iniloc <- data.frame(matrix(c(lubridate::day(meta$time_coverage_start[i]),
                                lubridate::month(meta$time_coverage_start[i]),
                                lubridate::year(meta$time_coverage_start[i]),
                                meta$geospatial_lat_start[i],
                                make360(meta$geospatial_lon_start[i]),
                                lubridate::day(meta$time_coverage_end[i]),
                                lubridate::month(meta$time_coverage_end[i]),
                                lubridate::year(meta$time_coverage_end[i]),
                                meta$geospatial_lat_end[i],
                                make360(meta$geospatial_lon_end[i])),
                              nrow = 2, ncol = 5, byrow = T))
  
  names(iniloc) <- list('day','month','year','lat','lon')
  iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
  tag <- iniloc$date[1]
  pop <- iniloc$date[2]
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  
  sp.lim <- list(lonmin = make360(limits$lonmin[which(limits$instrument_name == meta$instrument_name[i])]),
                 lonmax = make360(limits$lonmax[which(limits$instrument_name == meta$instrument_name[i])]),
                 latmin = limits$latmin[which(limits$instrument_name == meta$instrument_name[i])],
                 latmax = limits$latmax[which(limits$instrument_name == meta$instrument_name[i])])
  
  ## setup the spatial grid to base likelihoods on
  locs.grid <- setup.locs.grid(sp.lim, res='quarter')
  
  work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
  setwd(work_dir)
  
  ## cut the sst and ohc data from a subset of individual's likelihoods due to external thermistor failure
  if (length(grep('B2381', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2005-05-11', tz='UTC')
  if (length(grep('D1464', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2007-04-29', tz='UTC')
  if (length(grep('A1973', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2004-08-14', tz='UTC')
  if (length(grep('A2088', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2004-12-01', tz='UTC')
  
  if (length(grep('390191', meta$instrument_name[i])) > 0){
    cut_date <- as.POSIXct('2004-01-19', tz='UTC')
    iniloc$day[2] <- 19
    iniloc$month[2] <- 01
    iniloc$lat[2] <- 24.248
    iniloc$lon[2] <- 243.64
    iniloc$date[2] <- cut_date
    pop <- cut_date
    dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  }
  
  data_dir <- paste0('~/ebs/Data/albacore/', meta$instrument_name[i], '/cdb/')
  #  data_dir <- paste0('~/Google Drive File Stream/My Drive/Albacore - All Data/data/', meta$instrument_name[i], '/cdb/')
  etuff_file <- paste(data_dir, meta$instrument_name[i], '_eTUFF.txt', sep='')
  
  etuff <- read_archival(etuff_file)
  
  ## sst is Date, Temperature
  ## most tags have some reported SST metric but easier to create our own with archival tags
  #sstVars <- varNames[grep('sst', varNames)]
  #sst1 <- archival_to_etuff(etuff$etuff, vars = c('DateTime', sstVars))
  
  sst2 <- archival_to_etuff(etuff$etuff, vars ='sst')
  rm(etuff)
  sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T), mean_sst = mean(sstMean, na.rm=T))
  #sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T))
  sst2$Date <- as.POSIXct(sst2$`as.Date(DateTime)`, tz='UTC')
  sst2 <- data.frame(sst2)
  for (b in 1:nrow(sst2)){
    if(is.nan(sst2$Temperature[b]) & !is.nan(sst2$mean_sst[b])) sst2$Temperature[b] <- sst2$mean_sst[b]
  }
  sst2 <- sst2 %>% dplyr::select(Date, Temperature)
  sst.dir <- paste0('./tmp/sst/')
  L.sst <- calc.sst.par(sst2, filename='oi', sst.dir = sst.dir, dateVec = dateVec, sens.err = 3, focalDim = 3)
  writeRaster(L.sst, filename = paste0('L.sst_', format(Sys.Date(), '%Y%m%d'), '.grd'), overwrite=TRUE)
  L.sst <- brick(paste0('L.sst_', format(Sys.Date(), '%Y%m%d'), '.grd'))
  #L.sst <- readAll(L.sst)
  
  fList <- list.files()
  load(file=fList[max(grep('L.res', fList))])
  
  L.res$L.rasters$L.sst <- readAll(raster::resample(L.sst, L.res$L.rasters$L.sst))
  save(L.res, file=paste0(meta$instrument_name[i],'_L.res_', format(Sys.Date(), '%Y%m%d'), '.rda'))
  rm(L.sst); gc()
  
  #if (i == 1 | i == 2) L.res$L.rasters$L.ohc <- flip(L.res$L.rasters$L.ohc, 'y')
=======
  # RUN THE FILTER STEP
  K <- list(K1, K2)
  f <- hmm.filter(g = L.res$g, L = L, K = K, maskL = FALSE, P = P, m = 2)
  nllf <- -sum(log(f$psi[f$psi > 0])) # negative log-likelihood
  aic <- 2 * nllf + 2 * length(which(!is.na(pars)))
  
  # RUN THE SMOOTHING STEP
  s <- hmm.smoother(f, K = K, L = L, P = P)
  
  # GET THE MOST PROBABLE TRACK
  tr <- calc.track(s, g = L.res$g, dateVec, iniloc)
  plotHMM(s, tr, dateVec, ptt=paste(tt), save.plot = F)
  
  # WRITE OUT RESULTS
  outVec <- matrix(c(tt, pars, aic = aic), ncol=6)
  res <- list(outVec = outVec, s = s, g = L.res$g, tr = tr, dateVec = dateVec, iniloc = iniloc)
  save(res, file=paste0(meta$instrument_name[i], '_', tt, '_HMMoce_res_', format(Sys.Date(), '%Y%m%d'), '.rda'))
>>>>>>> 2597194652866eb91bd5bf3a16df89994df86710
  
  write.table(res$outVec, file=paste0(meta$instrument_name[i], '_HMMoce_results_', format(Sys.Date(), '%Y%m%d'), '.csv'), 
              sep=',', append=T, row.names = F, col.names = F)
  
<<<<<<< HEAD
  if (exists('cut_date')){
    for (b in which(dateVec >= cut_date)) L.res$L.rasters$L.sst[[b]] <- L.res$L.rasters$L.sst[[b]] * NA
    for (b in which(dateVec >= cut_date)) L.res$L.rasters$L.ohc[[b]] <- L.res$L.rasters$L.ohc[[b]] * NA
  }
   
  ## each run idx combination
  for (tt in 1:length(run_idx)){
    
    #iniloc$lon <- make360(iniloc$lon)
    t1 <- Sys.time()
    L <- HMMoce:::make.L.par(L.res$L.rasters[run_idx[[tt]]], iniloc, dateVec, 
                             ncores = round(parallel::detectCores() / 2, 2))
    t2 <- Sys.time()
    print(t2-t1)
    print(tt)
    #[1] "Processing in parallel using 36 cores... " m4.10xlarge
    # Time difference of 12.06928 mins for 0396
    if (length(grep('390191', meta$instrument_name[i])) > 0) L <- L[1:177,,]
    save(L, file=paste0(meta$instrument_name[i],'_L_tt', tt, '_', format(Sys.Date(), '%Y%m%d'), '.rda'))
    
    #load(fList[grep('L_tt3', fList)])
    gc()
    
    #L.mle <- coarse.L(L, L.res$L.rasters)$L.mle
    #g.mle <- coarse.L(L, L.res$L.rasters)$g.mle
    #pars.optim <- HMMoce:::opt.params(pars.init = c(2,.2,.6,.8), 
    #                         lower.bounds = c(0.1, 0.001, .1, .1),  
    #                         upper.bounds = c(6, .6, .9, .9), 
    #                         g = L.res$g, 
    #                         #g = g.mle, 
    #                         L = L, 
    #                         #L = L.mle, 
    #                         alg.opt = 'optim', 
    #                         write.results = FALSE)
    #pars.optim.mle <- opt.params(pars.init = c(2,.2,.6,.8), 
    #                             lower.bounds = c(0.1, 0.001, .1, .1),  
    #                             upper.bounds = c(6, .6, .9, .9), 
    #                             #g = L.res$g, 
    #                             g = g.mle, 
    #                             #L = L, 
    #                             L = L.mle, 
    #                             alg.opt = 'optim', 
    #                             write.results = FALSE)
    
    ## about 22 mins per iteration on blue shark 141259
    ## about 1.5 mins with MLE grid
    ## final vals way different between the different grids
    
    #pars.nlminb <- opt.params(pars.init = c(2,.2,.6,.8), 
    #                          lower.bounds = c(0.1, 0.001, .1, .1),  
    #                          upper.bounds = c(6, .6, .9, .9), 
    #                          g = L.res$g, 
    #                          #g = g.mle, 
    #                          L = L, 
    #                          #L = L.mle, 
    #                          alg.opt = 'nlminb', 
    #                          write.results = FALSE)
    ## about 30 mins per iteration on blue shark 141259
     
    #ceiling(parallel::detectCores() * .9)
    #pars.ga.one <- opt.params(pars.init = c(2), 
    #                          lower.bounds = c(1), 
    #                          upper.bounds = c(8), 
    #                          g = L.res$g, 
    #                          #g = g.mle, 
    #                          L = L, 
    #                          #L = L.mle, 
    #                          alg.opt = 'ga', 
    #                          write.results = FALSE,
    #                          ncores = ceiling(parallel::detectCores() * .9))
    #ncores = 2)
    
    #pars.ga.mle <- opt.params(pars.init = c(2,.2,.6,.8), 
    #                          lower.bounds = c(0.1, 0.001, .1, .1), 
    #                          upper.bounds = c(6, .6, .9, .9), 
    #                          #g = L.res$g, 
    #                          g = g.mle, 
    #                          #L = L, 
    #                          L = L.mle, 
    #                          alg.opt = 'ga', 
    #                          write.results = FALSE,
    #                          ncores = ceiling(parallel::detectCores() * .9))
    #ncores = 2)
    
    pars.ga <- HMMoce:::opt.params(pars.init = c(2,.2,.6,.8), 
                          lower.bounds = c(0.1, 0.001, .1, .1), 
                          upper.bounds = c(6, .6, .9, .9),
                          max_iter = 100,
                          run = 10,
                          g = L.res$g, 
                          #g = g.mle, 
                          L = L, 
                          #L = L.mle, 
                          alg.opt = 'ga', 
                          write.results = FALSE,
                          ncores = 10)
                          #ncores = ceiling(parallel::detectCores() / 2))
    #ncores = ceiling(parallel::detectCores() * .9))
    # 10 out of 40 cores on m4.10xl
    # Time difference of 7.049704 hours for 0396
    # [1,] 5.039326 0.1911284 0.7746087 0.8083944
    #Time difference of 1.393883 days for 1090269
    
    pars <- pars.ga$par
    sigmas = pars[1:2]
    sizes = rep(ceiling(sigmas[1]*4),2)
    pb = pars[3:4]
    muadvs = c(0,0)
    
    ## behav 1
    if(sizes[1]%%2==0){sizes[1]=sizes[1]+1}
    ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
    while(ss<.999){
      sizes[1]=sizes[1]+2
      ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
    }
    K1 <- HMMoce:::gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
    rm(ss)
    K1 <- HMMoce:::mask.K(K1)
    
    ## behav 2
    if(sizes[2] %% 2 == 0) sizes[2] <- sizes[2] + 1
    ss <- sum(gausskern.nostd(sizes[2], sigmas[2], muadv = muadvs[2]))
    while(ss < .999){
      sizes[2] <- sizes[2] + 2
      ss <- sum(gausskern.nostd(sizes[2], sigmas[2], muadv = muadvs[2]))
    }
    K2 <- HMMoce:::gausskern.pg(sizes[2], sigmas[2], muadv = muadvs[2])
    K2 <- HMMoce:::mask.K(K2)
    
    P <- matrix(c(pb[1], 1 - pb[1], 1 - pb[2], pb[2]), 2, 2, byrow = TRUE)
    
    # RUN THE FILTER STEP
    K <- list(K1, K2)
    f <- hmm.filter(g = L.res$g, L = L, K = K, maskL = FALSE, P = P, m = 2)
    nllf <- -sum(log(f$psi[f$psi > 0])) # negative log-likelihood
    aic <- 2 * nllf + 2 * length(which(!is.na(pars)))
     
    # RUN THE SMOOTHING STEP
    s <- hmm.smoother(f, K = K, L = L, P = P)
    
    # GET THE MOST PROBABLE TRACK
    tr <- calc.track(s, g = L.res$g, dateVec, iniloc)
    plotHMM(s, tr, dateVec, ptt=paste(tt), save.plot = F)
    
    # WRITE OUT RESULTS
    outVec <- matrix(c(tt, pars, aic = aic), ncol=6)
    res <- list(outVec = outVec, s = s, g = L.res$g, tr = tr, dateVec = dateVec, iniloc = iniloc)
    save(res, file=paste0(meta$instrument_name[i], '_', tt, '_HMMoce_res_', format(Sys.Date(), '%Y%m%d'), '.rda'))
    
    write.table(res$outVec, file=paste0(meta$instrument_name[i], '_HMMoce_results_', format(Sys.Date(), '%Y%m%d'), '.csv'), 
                sep=',', append=T, row.names = F, col.names = F)
    
     
    rm(s); rm(f); rm(pars.ga); rm(pars); rm(tr); rm(outVec); rm(L)
    gc()
    #load(file=past e0(meta$instrument_name[i],'_L.res_', format(Sys.Date(), '%Y%m%d'), '.rda'))
    
     
  } ## run_idxa
  
  rm(L.rasters); rm(L.res); rm(L); rm(cut_date); gc()
  #file.remove(list.files(hycom.dir, full.names = T))
  #file.remove(list.files(sst.dir, full.names = T))
  
} ## individual loop
=======
} ## run_idx

file.remove(list.files(hycom.dir, full.names = T))
>>>>>>> 2597194652866eb91bd5bf3a16df89994df86710

} ## individual loop

## this should be the end of the model building and we're ready for model selection?

## this should be the end of the model building and we're ready for model selection?




for (i in 1:nrow(meta)){
  
  work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
  setwd(work_dir)
  
  fList <- list.files()
  res_files <- fList[grep('results', fList)]
  res_files <- res_files[grep('.csv', res_files)]
  
  for (b in 1:length(res_files)){
    
    df <- read.table(res_files[b], sep=',', header=F)
    
    ## get run date
    df$run_date <- substr(res_files[b], 
                          stringr::str_locate(res_files[b], 'results_')[1,2] + 1,
                          nchar(res_files[b]) - 4)
    
    df$instrument_name <- meta$instrument_name[i]
    
    write.table(df, file='../albacore_hmmoce_runs.csv', append=TRUE, sep=',', col.names = F, row.names = F)
    
  }
  
  if (!file.exists(paste0(base_dir,'/sync_albacore_aws.sh'))){
    write.table('#!/bin/bash',
                file=paste0(base_dir,'/sync_albacore_aws.sh'), sep='', col.names = F, row.names = F, quote=F, append=F)
  }
  write.table(paste0('aws s3 sync ./', meta$instrument_name[i], '/ s3://braunmpg/RCode/HMMoce_run/', meta$instrument_name[i], '/'),
              file=paste0(base_dir,'/sync_albacore_aws.sh'), sep='', col.names = F, row.names = F, quote=F, append=T)
  
  
}



raw <- raw %>% filter(date <= iniloc$date[2])

raw[nrow(raw) + 1,] <- c('160424_2015_141261.1', NA, '3', iniloc$lon[2], iniloc$lat[2])
raw$date[nrow(raw)] <- iniloc$date[2]
raw$lon <- as.numeric(raw$lon); raw$lat <- as.numeric(raw$lat)

df.locs <- split(raw, raw$id)
ssm_fit <- lapply(df.locs, function(x) foieGras::fit_ssm(x, model='crw', time.step = 24, vmax=10))

## grab predicted locations output from fit_ssm
plocs <- lapply(ssm_fit, FUN=function(x) foieGras::grab(x, what = "p", as_sf = FALSE)) %>%
  do.call(rbind, .) %>%
  tbl_df() %>%
  mutate(id = as.character(id)) %>% group_by(id)

## subsample from predicted locations, if applicable
res_out <- 24
if (!is.null(res_out)){
  plocs <- split(plocs, plocs$id)
  
  for (i in 1:length(plocs)){
    res.i <- as.numeric(difftime(plocs[[i]]$date[2], plocs[[i]]$date[1], units = 'hours'))
    res_slice <- res_out / res.i
    if (res_slice < 1) res_slice <- 1
    plocs[[i]] <- plocs[[i]] %>% do(slice(., seq(1, n(), by = res_slice)))
  }
  
  plocs <- plocs %>%
    do.call(rbind, .) %>%
    tbl_df() %>%
    mutate(id = as.character(id)) %>% group_by(id)
  
}

mpm_fit <- try(fit_mpm(plocs[,c('id','date','lon','lat')], model = 'mpm'), silent=TRUE)
plocs$g <- NA
for (tt in 1:length(mpm_fit$mpm)){
  
  if (mpm_fit$converged[tt]){
    temp_g <- mpm_fit$mpm[[tt]]$fitted$g
    plocs$g[which(plocs$id == unique(plocs$id)[tt])] <- temp_g
    
  } else {
    #temp_g <- rep(NA, length.out = nrow(plocs[which(plocs$id == unique(plocs$id)[tt]),]))
  }
  
}

## convert to same scale as hmmoce behav estimate "p"
plocs$p_equiv <- plocs$g * -1 + 1
plocs$dateVec <- findInterval(plocs$date, dateVec)


out <- read.table('141259_HMMoce_results_outVec_20200824.csv', sep=',', header=F)
names(out) <- c('likvec','pars','ii','pars1','pars2','pars3','pars4','nll')
out[c((ncol(out) + 1):(ncol(out) + 8))] <- NA
names(out)[9:ncol(out)] <- c('rmse.lon','rmse.lat','gcd_mean','gcd_sd','gcd_median','gcd_min','gcd_max','rmse.behav')
out$likvec <- as.character(out$likvec)

for (i in 1:nrow(out)){
  #stringr::str_locate(out$likvec[i], '_')
  if (nchar(out$likvec[i]) > 2){
    result <- try(load(paste0(substr(out$likvec[i],1,1), substr(out$likvec[i],3,5), '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200828.rda')), silent=TRUE)
    if (class(result) == 'try-error'){
      result <- try(load(paste0(substr(out$likvec[i],1,1), substr(out$likvec[i],3,5), '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200824.rda')), silent=TRUE)
      if (class(result) == 'try-error'){
        result <- try(load(paste0(substr(out$likvec[i],1,1), substr(out$likvec[i],3,5), '_', out$pars[i], '_', out$ii[i], '-HMMoce_res.rda')), silent=TRUE)
      }
    }
  } else{
    result <- try(load(paste0(out$likvec[i], '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200828.rda')), silent=TRUE)
    if (class(result) == 'try-error'){
      result <- try(load(paste0(out$likvec[i], '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200824.rda')), silent=TRUE)
      if (class(result) == 'try-error'){
        result <- try(load(paste0(out$likvec[i], '_', out$pars[i], '_', out$ii[i], '-HMMoce_res.rda')), silent=TRUE)
      }
    }
  } 
  
  if (class(result) == 'try-error') next
  # res <- list(outVec = outVec, s = s, g = L.res$g, tr = tr, dateVec = dateVec, iniloc = iniloc)
  
  comp <- compareTracks(res$tr, plocs, dateVec)
  
  out[i,9:ncol(out)] <- unlist(comp)
  
  rm(res)
}

summary(out)

out %>%
  select(-gcd_min, -gcd_max) %>%
  gather(-likvec, -pars, -ii, -pars1, -pars2, -pars3, -pars4, -nll, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = nll, color = factor(likvec), shape = factor(pars))) +
  geom_point() + ylim(90,200) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

out %>%
  select(pars, rmse.lon) %>%
  filter(pars %in% c(1,3,5)) %>%
  #gather(-likvec, -pars, -ii, -pars1, -pars2, -pars3, -pars4, -nll, key = "var", value = "value") %>% 
  ggplot(aes(as.factor(pars), rmse.lon)) +
  geom_boxplot() 

out %>%
  select(likvec, rmse.lon) %>%
  #gather(-likvec, -pars, -ii, -pars1, -pars2, -pars3, -pars4, -nll, key = "var", value = "value") %>% 
  ggplot(aes(as.factor(likvec), rmse.lon)) +
  geom_boxplot()

for (i in 1:length(dateVec)){
  image.plot(res$g$lon[1,], res$g$lat[,1], res$s[1,i,,])
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  invisible(readline(prompt=paste('Check plots. If passes QC, press [enter] to continue.', sep='')))
  
}




setwd('~/ebs/Data/HMMoce_run/141259/')
dir <- getwd()
#load('./141259_tryHMMoce_20200724.rda')
#save(L.res, file='141259_L.res_20200724.rda')
#save(plocs, file='141259_plocs_20200724.rda')

load('141259_L.res_20200724.rda')
load('141259_plocs_20200724.rda')

# set lik colors
lik.breaks <- seq(0, 1, length.out = 25)
lik.mid = lik.breaks[1:(length(lik.breaks)-1)]
#lik.col = jet.colors(length(lik.breaks)-1) #[as.vector((dataT))]
lik.col = terrain.colors(length(lik.breaks)-1, rev=T) #[as.vector((dataT))]

#png(paste0(dir, '/141259_diagnose.png'),
#    width=5000, height=3000, res=300, onefile=TRUE)
tt=12; bb=1; ii=1
load(paste(paste(tt, bb, ii, sep='_'), '-HMMoce_res.rda', sep=''))
iniloc <- res$iniloc; dateVec <- res$dateVec
tr <- res$tr

likVec <- c(1:5) 
if (length(likVec) > 2){
  combine_idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  combine_idx <- utils::combn(likVec, 2, simplify=F)
}


L <- make.L(L.res$L.rasters[combine_idx[[tt]]], iniloc, dateVec)

pars <- res$outVec[4:7]
sigmas=pars[1:2]
sizes=ceiling(sigmas*4)
pb=pars[3:4]
muadvs=c(0,0)

gausskern.PG <- ifelse(ii == 1, TRUE, FALSE)

if(gausskern.PG){
  # behav 1
  if(sizes[1]%%2==0){sizes[1]=sizes[1]+1}
  ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
  while(ss<.999){
    sizes[1]=sizes[1]+2
    ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
  }
  K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
  rm(ss)
  K1=mask.K(K1)
  
  # behav 2
  if(sizes[2]%%2==0){sizes[2]=sizes[2]+1}
  ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
  while(ss<.999){
    sizes[2]=sizes[2]+2
    ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
  }
  K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
  rm(ss)
  K2=mask.K(K2)
  
}else{
  sizes=rep(ceiling(sigmas[1]*4),2)
  K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
  K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
}


P <- matrix(c(pb[1],1-pb[1],1-pb[2],pb[2]),2,2,byrow=TRUE)


f <- hmm.filter(L.res$g, L, K1, K2, maskL=FALSE, P)

f.sum <- apply(f$pred, 2:4, sum)
s.sum <- apply(res$s, 2:4, sum)


pdf(paste0(dir, '/', paste0(paste(tt, bb, ii, sep='_'), '_diagnose.pdf')),
    width=12, height=8, onefile=TRUE)

#for (i in 1:10){
for (i in 1:length(dateVec)){
  
  if (length(combine_idx[[tt]]) == 2) nf <- layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), nrow=6, ncol=2, byrow=F), widths=c(5,5), heights=c(4,4,4,4,4,4))
  if (length(combine_idx[[tt]]) == 3) nf <- layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow=F), widths=c(5,5), heights=c(4,4,4))
  if (length(combine_idx[[tt]]) == 4) nf <- layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,6,7,7), nrow=6, ncol=3, byrow=F), widths=c(5,5,5), heights=c(4,4,4,4,4,4))
  
  #layout.show(nf)
  
  #========
  ## plot likelihood 1
  #========
  par (mar=c(2,4,4,2))
  
  image(L.res$L.rasters[combine_idx[[tt]]][[1]][[i]], col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[combine_idx[[tt]][[1]]])
  axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
  
  
  #========
  ## plot likelihood 2
  #========
  par (mar=c(2,4,4,2))
  
  image(L.res$L.rasters[combine_idx[[tt]]][[2]][[i]], col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes = F, main=names(L.res$L.rasters)[combine_idx[[tt]][[2]]])
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
  # legend
  #image(1, lik.mid, t(as.matrix(lik.mid)), breaks=lik.breaks, col=lik.col, axes=FALSE, xlab="",
  #      ylab=parse(text=paste('Temperature', "*degree~C", sep="")))
  #axis(2);box();
  
  if (length(combine_idx[[tt]]) >= 3){
    #========
    ## plot likelihood 3
    #========
    
    image(L.res$L.rasters[combine_idx[[tt]]][[3]][[i]], col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[combine_idx[[tt]][[3]]])
    axis(1); axis(2); box()
    world(add=T)
    points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    
  }
  
  
  
  if (length(combine_idx[[tt]]) == 4){
    #========
    ## plot likelihood 4
    #========
    
    image(L.res$L.rasters[combine_idx[[tt]]][[4]][[i]], col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[combine_idx[[tt]][[4]]])
    axis(1); axis(2); box()
    world(add=T)
    points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    
  }
  
  
  #========
  ## plot L
  #========
  image(L.res$g$lon[1,], L.res$g$lat[,1], L[i,,], col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes=F, main=dateVec[i])
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
  #========
  ## plot filter
  #========
  image(L.res$g$lon[1,], L.res$g$lat[,1], f.sum[i,,] / max(f.sum[i,,]), col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes=F, main='filter')
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
  
  
  #========
  ## plot smoother
  #========
  image(L.res$g$lon[1,], L.res$g$lat[,1], s.sum[i,,] / max(s.sum[i,,]), col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes=F, main='smoother')
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
}
dev.off()




#================
## FILTER LIGHT-BASED LON ESTIMATES
#================
# This is already done for WC tags as they've all been filtered via GPE2
# Need to filter the Lotek "raw" estimates

drop_light <- list()
for (i in 1:nrow(meta)){
  
  if (meta$manufacturer[i] != 'Lotek'){
    drop_light[[i]] <- NA
  } else{
    # light-based locations
    data_dir <- paste0('~/Google Drive File Stream/My Drive/Albacore - All Data/data/', meta$instrument_name[i], '/cdb/')
    etuff_file <- paste(data_dir, meta$instrument_name[i], '_eTUFF.txt', sep='')
    etuff <- read_archival(etuff_file)
    tr <- get_track(etuff)
    tr$lon2 <- make360(tr$longitude)
    
    ## do a little bit of automatic filtering and keep an index of those positions
    auto_idx <- which(tr$lon2 < 125 | tr$lon2 > 260)
    if (length(drop_auto) > 0){
      drop_auto <- tr[auto_idx,]
      tr <- tr[-auto_idx,]
    }
    
    world <- map_data('world2')
    xl <- c(min(tr$lon2) - 2, max(tr$lon2) + 2)
    yl <- c(min(tr$latitude) - 2, max(tr$latitude) + 2)
    
    ## simple map of move data
    m1 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
      coord_fixed(xlim=xl, ylim=yl) + xlab('') + ylab('') +
      geom_point(data = tr, aes(x = lon2, y = latitude, colour = DateTime))
    
    m2 <- ggplot() + geom_point(data = tr, aes(x = lon2, y = DateTime, colour=DateTime)) +
      geom_point(data = meta[i,], aes(x=make360(geospatial_lon_start), y=time_coverage_start), col='green') +
      geom_point(data = meta[i,], aes(x=make360(geospatial_lon_end), y=time_coverage_end), col='red')
    
    #m3 <- ggplot() + geom_point(data = tr, aes(x = DateTime, y = latitude, colour=DateTime))
    
    lay <- rbind(c(1,2))
    g <- gridExtra::marrangeGrob(grobs = list(m1, m2), heights = c(8),
                                 width = c(5), layout_matrix = lay)
    g
    
    plot(tr$lon2, tr$latitude, col='red'); 
    plot(tr$lon2, tr$DateTime, ylim=c(meta$time_coverage_start[i], meta$time_coverage_end[i])); 
    points(make360(meta$geospatial_lon_start[i]), meta$time_coverage_start[i], pch=24, bg='green')
    points(make360(meta$geospatial_lon_end[i]), meta$time_coverage_end[i], pch=23, bg='red')
    #fields::world(add=T, wrap=c(0,360))#, xlim=c(min(tr$lon2), max(tr$lon2)), ylim=c(min(tr$latitude), max(tr$latitude)))
    by_hand <- identify(tr$lon2, tr$DateTime)
    
    drop_list[[i]] <- list(drop_auto = drop_auto, by_hand = by_hand)
  }
  
}


<<<<<<< HEAD

## assigning pseudo end location for 390191
## visually inspected:

par(mfrow=c(2,2))
plot(L.light[[176]])
plot(L.sst[[176]])
plot(L.ohc[[176]])
plot(L.bathy[[176]])

## then made the combined likelihood L and visualized to find likely end location
L <- make.L(xxx)
image.plot(L[176,,]) ## annoying array format...

crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
lon.agg <- seq(min(locs.grid$lon[1,]), max(locs.grid$lon[1,]), 
               length.out = length(locs.grid$lon[1,]))
lat.agg <- seq(min(locs.grid$lat[,1]), max(locs.grid$lat[,1]), 
               length.out = length(locs.grid$lat[,1]))
list.L <- list(x = lon.agg, y = lat.agg, z = L[176,,])
ex <- raster::extent(list.L)
r <- raster::raster(list.L$z, xmn = ex[1], xmx = ex[2], 
                    ymn = ex[3], ymx = ex[4], crs)
r <- raster::flip(r, direction = "y")
plot(r)

## lon
.73*(248.5-230.5) + 230.5
# [1] 243.64

## lat 
.32 * (42.2 - 15.8) + 15.8
# [1] 24.248


#========================
#========================

## aggregate individual run results spreadsheets into master


for (i in 1:nrow(meta)){
  
  work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
  setwd(work_dir)
  
  fList <- list.files()
  results_idx <- grep('HMMoce_results', fList)
  
  for (bb in results_idx){
    
    df <- read.table(fList[bb], sep=',', header=F)
    
    grab_date <- stringr::str_locate_all(fList[bb], 'results')
    grab_date <- substr(fList[bb], grab_date[[1]][2] + 2, nchar(fList[bb]) - 4)
    
    df$build_date <- grab_date
    df$id <- meta$instrument_name[i]
    
    write.table(df, file='~/ebs/Data/albacore/albacore_HMMoce_runs.csv', sep=',', append = TRUE, col.names = FALSE, row.names = FALSE)
    
    rm(df); rm(grab_date)
    
  }
  
}

df <- read.table('~/ebs/Data/albacore/albacore_HMMoce_runs.csv', sep=',', header=F)
names(df) <- c('tt','par1','par2','par3','par4','aic','build_date','instrument_name')
 
## translate tt into actual runs
# tt = 1: light + sst + ohc + bathy
# tt = 2: light + sst + bathy
# tt = 3: light + ohc + bathy
df2 <- rbind(filter(df, tt == 1) %>% mutate(L = 'all'),
             filter(df, tt == 2) %>% mutate(L = 'drop_ohc'),
             filter(df, tt == 3) %>% mutate(L = 'drop_sst'))

df <- df %>% select(instrument_name, everything())

ggplot(df2, aes(fill=L, y=aic, x=instrument_name)) + 
  geom_bar(position="dodge", stat="identity")


=======
>>>>>>> 2597194652866eb91bd5bf3a16df89994df86710
