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
#library(GA) ## for genetic algorithm

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

#for (i in 1:nrow(meta)){

for (i in 3:nrow(meta)){
  
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
                                meta$geospatial_lon_start[i],
                                lubridate::day(meta$time_coverage_end[i]),
                                lubridate::month(meta$time_coverage_end[i]),
                                lubridate::year(meta$time_coverage_end[i]),
                                meta$geospatial_lat_end[i],
                                meta$geospatial_lon_end[i]), nrow = 2, ncol = 5, byrow = T))
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
  if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
  for (b in 1:length(dateVec)){
    if (!file.exists(paste0(sst.dir, '/', 'oi_', dateVec[b], '.nc'))) get.env(dateVec[b], filename='oi', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
  }
  
  hycom.dir <- paste0('./tmp/hycom/')
  if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
  
  ## check for hycom. if its not available pull it down from the s3 bucket
  #if (length(list.files(hycom.dir)) != length(dateVec)){
  #  system(paste0('aws s3 cp s3://braun-data/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/ /home/rstudio/ebs/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/'))
  #  if (length(list.files(hycom.dir)) != length(dateVec)) stop('Check that all daily hycom data is available for this individual run.')
  #}
  
  ## this is how we get hycom. it will skip days for which data already exists
  for (b in 1:length(dateVec)){
    if (!file.exists(paste0(hycom.dir, '/', 'hycom_', dateVec[b], '.nc'))) {
      try(
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
  
} ## run_idx

file.remove(list.files(hycom.dir, full.names = T))

} ## individual loop


## this should be the end of the model building and we're ready for model selection?









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


