#================
## INTRO STUFF
#================

#devtools::load_all('../HMMoce')
#devtools::install_github('camrinbraun/HMMoce', ref='dev')
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

base_dir <- '~/work/RCode/HMMoce_run'
setwd(base_dir)

meta <- read.table('../nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, stringsAsFactors = F)
meta <- meta %>% filter(person_owner == 'Austin Gallagher' & instrument_type == 'popup' & manufacturer == 'Lotek')
meta$time_coverage_start <- as.POSIXct(meta$time_coverage_start, tz = 'UTC')
meta$time_coverage_end <- as.POSIXct(meta$time_coverage_end, tz = 'UTC')

for (i in 1:nrow(meta)){
  
  ## temporal bounds
  iniloc <- data.frame(matrix(c(lubridate::day(meta$time_coverage_start[i]),
                                lubridate::month(meta$time_coverage_start[i]),
                                lubridate::year(meta$time_coverage_start[i]),
                                meta$geospatial_lat_start[i],
                                meta$geospatial_lon_start[i],
                                lubridate::day(meta$time_coverage_end[i]),
                                lubridate::month(meta$time_coverage_end[i]),
                                lubridate::year(meta$time_coverage_end[i]),
                                meta$geospatial_lat_end[i],
                                meta$geospatial_lon_end[i]),
                              nrow = 2, ncol = 5, byrow = T))
  
  names(iniloc) <- list('day','month','year','lat','lon')
  iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
  tag <- iniloc$date[1]
  pop <- iniloc$date[2]
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  
  if (meta$instrument_name[i] == "160189_2019_183613"){
    limits <- list(lonmin=-85, lonmax=-70, latmin=18, latmax=30)
  }
  
  sp.lim <- limits
  
  ## setup the spatial grid to base likelihoods on
  locs.grid <- setup.locs.grid(sp.lim, res='hycom')
  
  work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
  if (!dir.exists(work_dir)) dir.create(work_dir)
  setwd(work_dir)
  
  sst.dir <- paste0('./tmp/sst/')
  if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
  for (b in 1:length(dateVec)){
    if (!file.exists(paste0(sst.dir, '/', 'ghr_', dateVec[b], '.nc'))) get.env(dateVec[b], filename='ghr', type = 'sst', sst.type='mur', spatLim = sp.lim, save.dir = sst.dir)
  }
  
  hycom.dir <- paste0('./tmp/hycom/')
  if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
  ## this is how we get hycom. it will skip days for which data already exists
  for (b in 1:length(dateVec)){
    if (!file.exists(paste0(hycom.dir, '/', 'hycom_', dateVec[b], '.nc'))) {
      print(b)
      try(
        get.env(dateVec[b], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir), TRUE)
    }
  }
  
  
  bathy.dir <- paste0('./tmp/bathy/')
  #if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
  #if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
  #  bathy <- HMMoce::get.bath.data(sp.lim, save.dir = bathy.dir, res=1)
  #} else{ ## OR (once downloaded and reading the .nc later)
  #bathy <- raster::raster('~/ebs/EnvData/bathy/global_bathy_0.01.nc')
  #bathy <- raster::crop(bathy, raster::extent(unlist(sp.lim)))
  #  raster::writeRaster(bathy, paste0(bathy.dir, '/bathy.nc'), format='CDF', overwrite=TRUE, varname = 'topo')                
  #}
  #print(i); gc()
  bathy <- raster::raster(paste0(bathy.dir, '/bathy.nc'))
  
  data_dir <- paste0('~/Desktop/data_org/', meta$instrument_name[i],'/')
  etuff_file <- paste(data_dir, meta$instrument_name[i], '_eTUFF.txt', sep='')
  
  etuff <- read_archival(etuff_file)
  
  ## sst is Date, Temperature
  sst2 <- archival_to_etuff(etuff$etuff, vars ='sst')
  sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T))#, mean_sst = mean(sstMean, na.rm=T))
  #sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T))
  sst2$Date <- as.POSIXct(sst2$`as.Date(DateTime)`, tz='UTC')
  sst2 <- data.frame(sst2)
  #for (b in 1:nrow(sst2)){
  #  if(is.nan(sst2$Temperature[b]) & !is.nan(sst2$mean_sst[b])) sst2$Temperature[b] <- sst2$mean_sst[b]
  #}
  sst2 <- sst2 %>% dplyr::select(Date, Temperature)
  
  ## mmd is Date, MaxDepth
  mmd2 <- archival_to_etuff(etuff$etuff, vars ='depthMax')
  mmd2 <- mmd2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), depthMax = max(depthMax, na.rm=T))#, mean_sst = mean(sstMean, na.rm=T))
  #mmd2 <- mmd2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T))
  mmd2$Date <- as.POSIXct(mmd2$`as.Date(DateTime)`, tz='UTC')
  mmd2 <- mmd2 %>% dplyr::select(Date, depthMax) %>% as.data.frame()
  names(mmd2)[grep('depthMax', names(mmd2))] <- 'MaxDepth'
  
  ## depth-temp profiles: Date, Depth, MeanTemp
  series <- get_series(etuff)
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
  #tr$lon2 <- make360(tr$longitude)
  auto_drop <- which(tr$longitude > -50 | tr$longitude < -100 | tr$latitude > 30)
  if (length(auto_drop) > 0) tr <- tr[-auto_drop,]
  #if (all(!is.na(drop_vec))) tr <- tr[-drop_vec,]
  ## then we won't have error estimates. Musyl error estimates will be applied automatically in the likelihood calc
  tr$Offset <- 0
  tr$Offset.orientation <- 0
  names(tr)[grep('longitudeError', names(tr))] <- 'Error.Semi.minor.axis'
  names(tr)[grep('latitudeError', names(tr))] <- 'Error.Semi.major.axis'
  tr$Error.Semi.major.axis <- tr$Error.Semi.major.axis * 110 * 1000 ## converts degrees error to meters
  tr$Error.Semi.minor.axis <- tr$Error.Semi.minor.axis * 110 * 1000
  names(tr)[grep('latitude', names(tr))] <- 'Latitude'
  names(tr)[grep('longitude', names(tr))] <- 'Longitude'
  names(tr)[grep('DateTime', names(tr))] <- 'Date'
  tr <- tr %>% dplyr::select(Date, Longitude, Latitude, Error.Semi.minor.axis, Error.Semi.major.axis, Offset, Offset.orientation)
  
  world <- map_data('world')
  xl <- c(min(tr$Longitude) - 5, max(tr$Longitude) + 5)
  yl <- c(min(tr$Latitude) - 5, max(tr$Latitude) + 5)
  
  ## simple map of move data just to double-check before calculating light-based likelihoods
  p1 <- ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group)) +
    coord_fixed(xlim=xl, ylim=yl) + xlab('') + ylab('') +
    geom_point(data=tr, aes(x=Longitude, y=Latitude, colour=Date)) +
    geom_point(data = meta[i,], aes(x=geospatial_lon_start, y=geospatial_lat_start), col='green') +
    geom_point(data = meta[i,], aes(x=geospatial_lon_end, y=geospatial_lat_end), col='red') #+
  p1
  
  ## LIGHT
  L.light <- calc.lightloc(tr, locs.grid = locs.grid, dateVec = dateVec, errEll = TRUE)

  ## SST
  #L.sst <- calc.sst.par(sst, filename='oi', sst.dir = sst.dir, dateVec = dateVec[1:10], sens.err = 3, focalDim = 3)
  L.sst <- calc.sst.par(sst2, filename='oi', sst.dir = sst.dir, dateVec = dateVec, sens.err = 3, focalDim = 3)
 
  # OCEAN HEAT CONTENT (INTEGRATED PDTS)
  L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec[3], isotherm = '', use.se = F)
  #L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, 
                        isotherm = '', use.se = F)#, ncores = 30)
  #writeRaster(L.ohc, filename = paste0(work_dir, 'L.ohc.grd'), format='raster')

  # HYCOM PROFILE BASED LIKELIHOODS
  #L.hycom <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:3], use.se = T)
  #L.hycom.se <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:3], use.se = T)
  L.hycom.par <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = F)#, ncores=2)
  
  
  L.bathy <- calc.bathy(mmd2, bathy, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')
  
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
                  c(1,3,4),
                  c(1,4),
                  c(1),
                  c(1:2))
  
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
    pars.ga.one <- HMMoce:::opt.params(pars.init = c(2), 
                                       lower.bounds = c(1), 
                                       upper.bounds = c(8), 
                                       g = L.res$g, 
                                       #                          #g = g.mle, 
                                       L = L, 
                                       #                          #L = L.mle, 
                                       alg.opt = 'ga', 
                                       write.results = FALSE,
                                       ncores = ceiling(parallel::detectCores() * .9))
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
    
    
    #pars <- pars.ga$par
    pars <- pars.ga.one$par
    if (length(pars) == 4){
      sigmas = pars[1:2]
      sizes = rep(ceiling(sigmas[1]*4),2)
      pb = pars[3:4]
      muadvs = c(0,0)
    } else if (length(pars) == 1){
      sigmas = pars[1]
      sizes = rep(ceiling(sigmas[1]*4),2)
      pb = NULL
      muadvs = c(0)
    }
    
    # behav 1
    if(sizes[1]%%2==0){sizes[1]=sizes[1]+1}
    ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
    while(ss<.999){
      sizes[1]=sizes[1]+2
      ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
    }
    K1 = HMMoce:::gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
    rm(ss)
    K1=HMMoce:::mask.K(K1)
    
    # behav 2
    if (!is.null(pb)){
      if(sizes[2]%%2==0){sizes[2]=sizes[2]+1}
      ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
      while(ss<.999){
        sizes[2]=sizes[2]+2
        ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
      }
      K2 = HMMoce:::gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
      rm(ss)
      K2=HMMoce:::mask.K(K2)
    }
    
    
    
    ## set transition matrix, if applicable
    if (!is.null(pb)){
      P <- matrix(c(pb[1],1-pb[1],1-pb[2],pb[2]),2,2,byrow=TRUE)
    } else{
      P <- NULL
    }
    
    
    # RUN THE FILTER STEP
    if (!is.null(pb)){
      K=list(K1,K2)
      f <- hmm.filter(g=L.res$g, L=L, K=K, maskL=FALSE, P=P, m=2)
    } else{
      K=list(K1)
      f <- hmm.filter(g=L.res$g, L=L, K=K, maskL=FALSE, P=P, m=1)
    }
    nllf <- -sum(log(f$psi[f$psi > 0])) # negative log-likelihood
    aic <- 2 * nllf + 2 * length(which(!is.na(pars)))
    
    
    # RUN THE SMOOTHING STEP
    s <- hmm.smoother(f, K = K, L = L, P = P)
    
    # GET THE MOST PROBABLE TRACK
    tr <- calc.track(s, g = L.res$g, dateVec, iniloc)
    #plotHMM(s, tr, dateVec, ptt=paste(tt), save.plot = F)
    
    # WRITE OUT RESULTS
    outVec <- matrix(c(tt, pars, aic = aic), ncol=6)
    res <- list(outVec = outVec, s = s, g = L.res$g, tr = tr, dateVec = dateVec, iniloc = iniloc)
    save(res, file=paste0(meta$instrument_name[i], '_', tt, '_HMMoce_res_', format(Sys.Date(), '%Y%m%d'), '.rda'))
    
    write.table(res$outVec, file=paste0(meta$instrument_name[i], '_HMMoce_results_', format(Sys.Date(), '%Y%m%d'), '.csv'), 
                sep=',', append=T, row.names = F, col.names = F)
    
    
    rm(s); rm(f); rm(pars.ga); rm(pars); rm(tr); rm(outVec); rm(L)
    gc()
    #load(file=past e0(meta$instrument_name[i],'_L.res_', format(Sys.Date(), '%Y%m%d'), '.rda'))
    
    
  } ## run_idx
  
  rm(L.rasters); rm(L.res); rm(L); rm(cut_date); gc()
  #file.remove(list.files(hycom.dir, full.names = T))
  #file.remove(list.files(sst.dir, full.names = T))
  
} ## individual loop

file.remove(list.files(hycom.dir, full.names = T))



