#========================
## HMMoce run w/example data
#========================
# might be a good idea to install latest version of HMMoce
# install.packages('HMMoce')
#devtools::install_github("camrinbraun/HMMoce", ref='dev', dependencies=F)
library(HMMoce); 
#devtools::load_all('~/ebs/RCode/HMMoce')
library(raster)
library(tidyverse)
# get necessary funs
#install.packages(c('readxl','readr')) #lubridate, plyr, MASS
library(readr); library(readxl); library(lubridate)
#devtools::install_github('camrinbraun/HMMoce', dependencies = F)
#source('https://raw.githubusercontent.com/galuardi/analyzepsat/v4.0/R/MWTextract.r')

## load all tag metadata master sheet
meta <- read.table('~/ebs/RCode/nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, blank.lines.skip = F, skip=0)
meta <- meta %>% filter(platform == 'Lamna nasus' & person_owner == 'Gregory Skomal' & ptt != 67829 & end_details != 'DNR')
meta$time_coverage_start <- as.Date(meta$time_coverage_start, format='%m/%d/%y', tz='UTC')
meta$time_coverage_end <- as.Date(meta$time_coverage_end, format='%m/%d/%y', tz='UTC')
udates <- seq(min(meta$time_coverage_start, na.rm=T), max(meta$time_coverage_end, na.rm=T), by='day')

#------------ 
##  GET ENVIRONMENTAL DATA 
#------------ 

# DOWNLOAD SST DATA
#sst.dir <- paste(tempdir(), '/sst/', sep='')
sst.dir <- '~/ebs/EnvData/sst/porbeagle/'
existing_sst <- as.Date(substr(list.files(sst.dir),8,17))
#dir.create(sst.dir, recursive = TRUE)
sst.udates <- udates[which(!(udates %in% existing_sst))]
get.env(sst.udates, filename='mursst', type = 'sst', sst.type='mur', spatLim = sp.lim, save.dir = sst.dir)

# YOU NEED SOME REPRESENTATION OF ENVIRONMENTAL DEPTH-TEMPERATURE
# HYCOM DATA
#hycom.dir <- paste(tempdir(), '/hycom/', sep='')
hycom.dir <- '~/ebs/EnvData/hycom/porbeagle/'
existing_hyc <- as.Date(substr(list.files(hycom.dir),7,16))
#dir.create(hycom.dir, recursive = TRUE)
pdt.udates <- udates[which(!(udates %in% existing_hyc))]
get.env(pdt.udates, filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

# BATHYMETRY
#bathy.dir <- paste(tempdir(), '/bathy/', sep='')
#bathy.dir <- '~/ebs/EnvData/bathy/porbeagle/'
#dir.create(bathy.dir, recursive = TRUE)
#bathy <- get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir, res=.5)
#bathy <- raster(paste(bathy.dir, 'bathy.nc', sep=''))
#bathy <- raster('~/ebs/EnvData/bathy/global_bathy.grd')
#bathy <- crop(bathy, extent(sp.lim$lonmin + 360, sp.lim$lonmax + 360, sp.lim$latmin, sp.lim$latmax))
#options(warn=0)
#bathy <- rotate(bathy)
#writeRaster(bathy, '~/ebs/EnvData/bathy/porbeagle/bathy_67829.grd')
bathy <- raster('~/ebs/EnvData/bathy/porbeagle/bathy_67829.grd')


#------------
# LOAD THE TAG DATA
#------------
row_idx <- which(meta$ptt %in% c(67826))
#row_idx <- which(meta$ptt %in% c(67824, 67827,67831,67833))
# SET INITIAL LOCATIONS (TAG AND POP-UP)
## we currently use a unique instrument name for each deployment (taxonserialnumber_deploymentyear_tagPTT)
for (idx in row_idx){
  #for (idx in 1:nrow(meta)){
  #idx <- i
  
  #idx <- grep('159911', meta$instrument_name)#which(meta$instrument_name == '159911_2006_67829')
  tag <- as.POSIXct(meta$time_coverage_start[idx], tz='UTC')
  pop <- as.POSIXct(meta$time_coverage_end[idx], tz='UTC')
  
  iniloc <- data.frame(matrix(c(lubridate::day(tag), lubridate::month(tag), lubridate::year(tag),
                                meta$geospatial_lat_start[idx], meta$geospatial_lon_start[idx], 
                                lubridate::day(pop), lubridate::month(pop), lubridate::year(pop),
                                meta$geospatial_lat_end[idx], meta$geospatial_lon_end[idx]),
                              nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  #tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  #pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day')) 
  
  # READ IN DATA AS OUTPUT FROM WC PORTAL
  # SST DATA
  sstFile <- paste('~/ebs/Data/data_org/', meta$instrument_name[idx], '/', meta$ptt[idx], '-SST.csv', sep='')
  tag.sst <- read.wc(sstFile, type = 'sst', tag=tag, pop=pop, verbose=T, dateFormat = '%H:%M:%S %d-%b-%Y') 
  sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data
  
  # DEPTH-TEMPERATURE PROFILE DATA
  pdtFile <- paste('~/ebs/Data/data_org/', meta$instrument_name[idx], '/', meta$ptt[idx], '-PDTs.csv', sep='')
  pdt <- read.wc(pdtFile, type = 'pdt', tag=tag, pop=pop, verbose=T) 
  pdt.udates <- pdt$udates; pdt <- pdt$data
  
  # RAW LIGHT DATA
  #lightFile <- system.file("extdata", "141259-LightLoc.csv", package = "HMMoce")
  #light <- read.wc(ptt, lightFile, type = 'light', tag=tag, pop=pop); 
  #light.udates <- light$udates; light <- light$data
  
  # LIGHT BASED POSITIONS FROM GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
  locsFile <- paste('~/ebs/Data/data_org/', meta$instrument_name[idx], '/', meta$ptt[idx], '-Locations-GPE2.csv', sep='')
  locs <- read.table(locsFile, sep = ',', header = T, blank.lines.skip = F)
  locDates <- as.Date(as.POSIXct(locs$Date, format=findDateFormat(locs$Date)))
  
  # SET SPATIAL LIMITS
  # these are the lat/lon bounds of your study area (e.g. where you think the animal went)
  sp.lim <- list(lonmin = -75,
                 lonmax = -55,
                 latmin = 35,
                 latmax = 47)
  
  
  #------------
  # CALCULATE LIKELIHOODS
  #------------
  # .par functions are the same calculations as those lacking .par, except they have been parallelized to leverage multiple CPUs
  locs.grid <- setup.locs.grid(sp.lim)
  
  # vector indicating which likelihoods to run (e.g. 1=light, 2=sst, 5=hycom)
  # can be combined with if() statements around calc functions: if (any(likVec == 5) & !exists('L.5')){calc.hycom(...)}
  likVec <- c(1,2,3,5) 
  
  # LIGHT-BASED LIKELIHOODS
  #L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25) # if trying to use raw light levels, not currently recommended (v0.2)
  L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
  #library(fields);library(raster)
  #plot(L.1[[12]]); world(add=T)
  
  # SST LIKELIHOODS
  L.2 <- calc.sst(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
  #L.2 <- calc.sst.par(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
  writeRaster(L.2, file=paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l2.grd', sep=''), overwrite=TRUE) # good idea to save after these larger calculations in case the next one causes problems
  #saveRDS(L.2, file=paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l2.RDS', sep='')) # good idea to save after these larger calculations in case the next one causes problems
  #L.2 <- readRDS(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l2.RDS', sep=''))
  gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
  
  # PDT LIKELIHOODS
  # OCEAN HEAT CONTENT (INTEGRATED PDTS)
  L.3 <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F, bathy = F)
  writeRaster(L.3, file=paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l3.grd', sep=''), overwrite=TRUE) # good idea to save after these larger calculations in case the next one causes problems
  #saveRDS(L.3, file=paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l3.RDS', sep='')) # good idea to save after these larger calculations in case the next one causes problems
  #L.3 <- readRDS(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l3.RDS', sep=''))
  gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
  
  # WORLD OCEAN ATLAS-BASED LIKELIHOODS
  #L.4 <- calc.woa.par(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
  # save.image() # good idea to save after these larger calculations in case the next one causes problems
  #gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
  
  # HYCOM PROFILE BASED LIKELIHOODS
  #pdt$Date <- as.POSIXct(pdt$Date, tz='UTC')
  L.5 <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
  writeRaster(L.5, file=paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l5.grd', sep=''), overwrite=TRUE) # good idea to save after these larger calculations in case the next one causes problems
  #saveRDS(L.5, file=paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l5.RDS', sep='')) # good idea to save after these larger calculations in case the next one causes problems
  #L.5 <- readRDS(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l5.RDS', sep=''))
  ## took 2 hrs on t2.large
  rm(L.2); rm(L.3); rm(L.5)
  gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
  
}

for (idx in row_idx[4]){
  #for (idx in 15:nrow(meta)){
  
  #idx <- grep('159911', meta$instrument_name)#which(meta$instrument_name == '159911_2006_67829')
  tag <- as.POSIXct(meta$time_coverage_start[idx], tz='UTC')
  pop <- as.POSIXct(meta$time_coverage_end[idx], tz='UTC')
  
  iniloc <- data.frame(matrix(c(lubridate::day(tag), lubridate::month(tag), lubridate::year(tag),
                                meta$geospatial_lat_start[idx], meta$geospatial_lon_start[idx], 
                                lubridate::day(pop), lubridate::month(pop), lubridate::year(pop),
                                meta$geospatial_lat_end[idx], meta$geospatial_lon_end[idx]),
                              nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  #tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  #pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day')) 
  
  # READ IN DATA AS OUTPUT FROM WC PORTAL
  # SST DATA
  sstFile <- paste('~/ebs/Data/data_org/', meta$instrument_name[idx], '/', meta$ptt[idx], '-SST.csv', sep='')
  tag.sst <- read.wc(sstFile, type = 'sst', tag=tag, pop=pop, verbose=T, dateFormat = '%H:%M:%S %d-%b-%Y') 
  sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data
  
  # DEPTH-TEMPERATURE PROFILE DATA
  pdtFile <- paste('~/ebs/Data/data_org/', meta$instrument_name[idx], '/', meta$ptt[idx], '-PDTs.csv', sep='')
  pdt <- read.wc(pdtFile, type = 'pdt', tag=tag, pop=pop, verbose=T) 
  pdt.udates <- pdt$udates; pdt <- pdt$data
  
  # RAW LIGHT DATA
  #lightFile <- system.file("extdata", "141259-LightLoc.csv", package = "HMMoce")
  #light <- read.wc(ptt, lightFile, type = 'light', tag=tag, pop=pop); 
  #light.udates <- light$udates; light <- light$data
  
  # LIGHT BASED POSITIONS FROM GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
  locsFile <- paste('~/ebs/Data/data_org/', meta$instrument_name[idx], '/', meta$ptt[idx], '-Locations-GPE2.csv', sep='')
  locs <- read.table(locsFile, sep = ',', header = T, blank.lines.skip = F)
  locDates <- as.Date(as.POSIXct(locs$Date, format=findDateFormat(locs$Date)))
  
  # SET SPATIAL LIMITS
  # these are the lat/lon bounds of your study area (e.g. where you think the animal went)
  sp.lim <- list(lonmin = -75,
                 lonmax = -55,
                 latmin = 35,
                 latmax = 47)
  
  
  #------------
  # CALCULATE LIKELIHOODS
  #------------
  # .par functions are the same calculations as those lacking .par, except they have been parallelized to leverage multiple CPUs
  locs.grid <- setup.locs.grid(sp.lim)
  
  # vector indicating which likelihoods to run (e.g. 1=light, 2=sst, 5=hycom)
  # can be combined with if() statements around calc functions: if (any(likVec == 5) & !exists('L.5')){calc.hycom(...)}
  likVec <- c(1,2,3,5) 
  
  L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
  #L.2 <- readRDS(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l2.RDS', sep=''))
  #L.3 <- readRDS(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l3.RDS', sep=''))
  #L.5 <- readRDS(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l5.RDS', sep=''))
  L.2 <- raster::brick(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l2.grd', sep=''))
  L.3 <- raster::brick(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l3.grd', sep=''))
  L.5 <- raster::brick(paste('~/ebs/Data/Porbeagles/', meta$ptt[idx], '_l5.grd', sep=''))
  
  #save.image(file='~/work/Data/thresher/196385/196385_hmmoce_likelihoods.rda') # good idea to save after these larger calculations in case the next one causes problems
  #gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
  #save.image('~/ebs/Data/AndersonJ/161211_Likelihoods.rda')
  
  #------------
  # PREPARE TO RUN THE MODEL
  #------------
  
  #load(paste('~/ebs/Data/AndersonJ/161211_Likelihoods.rda',sep=''))
  #known$lat <- as.numeric(as.character(known$lat))
  #known$lon <- as.numeric(as.character(known$lon))
  
  L.4 <- L.2 * 0 # just fill the missing likelihoods (that we're currently skipping) with zeros
  L.rasters <- mget(c('L.1', 'L.2', 'L.3', 'L.4', 'L.5')) # use with caution as all workspace items containing 'L.' will be listed. We only want the likelihood outputs calculated above
  #resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  resamp.idx <- 2
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  
  # Figure out appropriate L combinations
  # use this if you have a vector (likVec) indicating which likelihoods you are calculating
  # for example, likVec <- c(1,2,5) for light, sst, and hycom likelihoods
  if (length(likVec) > 2){
    L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
  } else{
    L.idx <- utils::combn(likVec, 2, simplify=F)
  }
  
  # which of L.idx combinations do you want to run?
  run.idx <- c(1,2,4,7,8)
  
  # vector of appropriate bounding in filter. see ?hmm.filter for more info
  bndVec <- c(NA, 1, 2)
  
  # vector of appropriate migr kernel speed. see ?makePar for more info.
  parVec <- c(1, 2)
  
  # GOOD IDEA TO CLEAN THINGS UP AND SAVE
  #rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))
  # setwd(); base::save.image('.rda')
  
  #------------
  # RUN THE MODEL
  #------------
  # CAN BE PARALLELIZED...
  #require(foreach)
  #print('Processing in parallel... ')
  #ncores <- ceiling(parallel::detectCores() * .25)
  #cl = parallel::makeCluster(ncores)
  #doParallel::registerDoParallel(cl, cores = ncores)
  #ans = foreach::foreach(tt = run.idx) %dopar%{
  ptt <- meta$ptt[idx]
  setwd(paste('~/ebs/Results/Porbeagles/',sep=''))
  pdt$dt <- pdt$Date
  pdt$Date <- as.Date(pdt$Date)
  for (tt in run.idx){
    for (bnd in bndVec){
      for (i in parVec){
        
        #ptt=141259
        runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i,sep='')
        
        # COMBINE LIKELIHOOD MATRICES
        # L.idx combination indicates likelihood surfaces to consider
        L <- make.L(L1 = L.res[[1]][L.idx[[tt]]],
                    L.mle.res = L.res$L.mle.res, dateVec = dateVec,
                    locs.grid = locs.grid, iniloc = iniloc, bathy = bathy,
                    pdt = pdt)#, known.locs = known)
        L.mle <- L$L.mle
        L <- L$L
        g <- L.res$g
        g.mle <- L.res$g.mle
        lon <- g$lon[1,]
        lat <- g$lat[,1]
        
        # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
        par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
        P.final <- par0$P.final
        
        # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
        par0 <- makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
        K1 <- par0$K1; K2 <- par0$K2
        
        # RUN THE FILTER STEP
        if(!is.na(bnd)){
          f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
          maskL.logical <- TRUE
        } else{
          f <- hmm.filter(g, L, K1, K2, P.final, maskL=F)
          maskL.logical <- FALSE
        }
        nllf <- -sum(log(f$psi[f$psi>0])) # negative log-likelihood
        
        # RUN THE SMOOTHING STEP
        s <- hmm.smoother(f, K1, K2, L, P.final)
        
        # GET THE MOST PROBABLE TRACK
        tr <- calc.track(s, g, dateVec, iniloc)
        #setwd(myDir); 
        plotHMM(s, tr, dateVec, ptt=runName, save.plot = T)
        
        # WRITE OUT RESULTS
        outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
                           Lidx = paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                           spLims = sp.lim[1:4], resol = raster::res(L.rasters[[resamp.idx]]),
                           maskL = maskL.logical, NLL = nllf, name = runName), ncol=15)
        #write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
        #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
        res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
        #setwd()
        save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        #save.image(file=paste(ptt, '-HMMoce.RData', sep=''))
        #source('~/HMMoce/R/hmm.diagnose.r') # not yet functional
        #hmm.diagnose(res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt, plot=T)
        
        write.table(outVec, file=paste(as.character(ptt), '_HMMoce_results_outVec.csv', sep=''), sep=',', append=T, row.names = F, col.names = F)
        
      } # parVec loop
    } # bndVec loop
  } # L.idx loop
  rm(L.rasters); rm(L.1); rm(L.2); rm(L.3); rm(L.4); rm(L.5); rm(L); removeTmpFiles(); gc()
}

#} # loop across indiv ptts
#parallel::stopCluster(cl)
#closeAllConnections()

fList <- list.files()
for (i in 1:nrow(meta)){
  
  print(as.character(meta$ptt[i]))
  print(length(fList[grep(meta$ptt[i], fList)]))
  
}



df <- read.table(paste(as.character(ptt), '_HMMoce_results_outVec.csv', sep=''), sep=',', header=F)

load('67829_idx2_bndNA_par2-HMMoce_res.rda') # load res

# run getCtr
bnds <- getCtr(res$s, res$tr, res$g, threshold = 5, makePlot=F)

# compile df of position, date, bnds (xtracto)
#data.frame(lapply(bnds, FUN=function(x) c(x$yDist, x$xDist)))
df <- cbind(res$tr, t(data.frame(lapply(bnds, FUN=function(x) c(x$yDist, x$xDist)))))
names(df)[5:6] <- c('ydist','xdist')
df$ydist[which(is.na(df$ydist))] <- mean(df$ydist, na.rm=T)
df$xdist[which(is.na(df$xdist))] <- mean(df$xdist, na.rm=T)


#ci <- list()
ci <- CI2shp(df)
plot(df$lon, df$lat, xlim=c(-72,-60), ylim=c(33,46), xlab='Lon', ylab='Lat')
plot(ci, add=T, col=alpha('grey', 0.5), border=NA)
points(df$lon, df$lat)
points(iniloc$lon[1], iniloc$lat[1], pch=16, col='green')
points(iniloc$lon[2], iniloc$lat[2], pch=16, col='red')

plot(rasterToContour(bathy, levels=-200), lty=2, add=T)

world(add=T, fill=T)



