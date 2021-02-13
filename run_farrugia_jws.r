#========================
## HMMoce run w/example data
#========================
# might be a good idea to install latest version of HMMoce
# install.packages('HMMoce')
devtools::install_github("camrinbraun/HMMoce", ref='dev', dependencies=F)
library(HMMoce); 
#devtools::load_all('~/work/RCode/HMMoce')
library(raster)
library(tidyverse)

# get necessary funs
#install.packages(c('readxl','readr')) #lubridate, plyr, MASS
#library(readr); library(readxl); library(lubridate)
#devtools::install_github('camrinbraun/HMMoce', dependencies = F)
#source('https://raw.githubusercontent.com/galuardi/analyzepsat/v4.0/R/MWTextract.r')

## load all tag metadata master sheet
meta <- read.table('~/Google Drive File Stream/My Drive/farrugia/add_etuff_JWS_SPOT_PAT_metadata.csv', sep=',', header=T, blank.lines.skip = F, skip=0)
meta <- meta %>% filter(instrument_type == 'popup')
meta$time_coverage_start <- as.POSIXct(meta$time_coverage_start, format='%m/%d/%y', tz='UTC')
meta$time_coverage_end <- as.POSIXct(meta$time_coverage_end, format='%m/%d/%y', tz='UTC')
ggplot(meta) + geom_segment(aes(x=time_coverage_start, xend=time_coverage_end, y=as.factor(ptt), yend=as.factor(ptt)))

#------------
# LOAD THE TAG DATA
#------------

for (idx in 1:nrow(meta)){
  
  work_dir <- paste0('~/work/RCode/HMMoce_run/', meta$ptt[idx], '/')
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
  setwd(work_dir)
  
  file.copy(paste0('~/Google Drive File Stream/My Drive/farrugia/', meta$ptt[idx], '/'), '~/work/RCode/HMMoce_run/', recursive = TRUE)
  
  # SET INITIAL LOCATIONS (TAG AND POP-UP)
  ## we currently use a unique instrument name for each deployment (taxonserialnumber_deploymentyear_tagPTT)
  tag <- as.POSIXct(meta$time_coverage_start[idx], format='%Y-%m-%d %H:%M:%S', tz='UTC')
  pop <- as.POSIXct(meta$time_coverage_end[idx], format='%Y-%m-%d %H:%M:%S', tz='UTC')
  
  iniloc <- data.frame(matrix(c(lubridate::day(tag), lubridate::month(tag), lubridate::year(tag),
                                meta$geospatial_lat_start[idx], meta$geospatial_lon_start[idx], 
                                lubridate::day(pop), lubridate::month(pop), lubridate::year(pop),
                                meta$geospatial_lat_end[idx], meta$geospatial_lon_end[idx]),
                              nrow = 2, ncol = 5, byrow = T))
  
  names(iniloc) <- list('day','month','year','lat','lon')
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day')) 
  
  # READ IN DATA AS OUTPUT FROM WC PORTAL
  fList <- list.files()
  
  # SST DATA
  sstFile <- fList[grep('-SST', fList)]
  tag.sst <- read.wc(sstFile, type = 'sst', tag=tag, pop=pop, verbose=T, dateFormat = '%H:%M:%S %d-%b-%Y') 
  sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data
  
  # DEPTH-TEMPERATURE PROFILE DATA
  pdtFile <- fList[grep('-PDT', fList)]
  pdt <- read.wc(pdtFile, type = 'pdt', tag=tag, pop=pop, verbose=T) 
  pdt.udates <- pdt$udates; pdt <- pdt$data
  
  # MAX DEPTH DATA
  mmdFile <- fList[grep('-MinMax', fList)]
  mmd <- read.table(mmdFile, sep=',', header=T)
  mmd$Date <- as.Date(as.POSIXct(mmd$Date, format='%H:%M:%S %d-%b-%Y', tz='UTC'))
  mmd <- mmd[which(!is.na(mmd$MaxDepth)),]
  mmd <- mmd[which(!duplicated(mmd$Date)),]
  
  # RAW LIGHT DATA
  #lightFile <- system.file("extdata", "141259-LightLoc.csv", package = "HMMoce")
  #light <- read.wc(ptt, lightFile, type = 'light', tag=tag, pop=pop); 
  #light.udates <- light$udates; light <- light$data
  
  # LIGHT BASED POSITIONS FROM GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
  locsFile <- fList[grep('-GPE2', fList)]
  locs <- read.table(locsFile, sep = ',', header = T, blank.lines.skip = F)
  locDates <- as.Date(as.POSIXct(locs$Date, format=findDateFormat(locs$Date)))
  
  # SET SPATIAL LIMITS
  # these are the lat/lon bounds of your study area (e.g. where you think the animal went)
  sp.lim <- rgdal::readOGR(paste0(meta$ptt[idx], '.shp'))
  sp.lim <- list(lonmin = extent(sp.lim)@xmin,
                 lonmax = extent(sp.lim)@xmax,
                 latmin = extent(sp.lim)@ymin,
                 latmax = extent(sp.lim)@ymax)
  
  #------------ 
  ##  GET ENVIRONMENTAL DATA 
  #------------ 
  # env data downloads can be
  #large, depending on application for 180 days of data spanning the NW Atlantic
  #(the example application), the downloads will take ~10mins on Amazon EC2.
  #Personal computers will likely be slower.
  
  # DOWNLOAD SST DATA
  #sst.dir <- paste(tempdir(), '/sst/', sep='')
  sst.dir <- './tmp/sst/'
  if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
  sst.udates <- unique(as.Date(tag.sst$Date))
  get.env(sst.udates, filename='mursst', type = 'sst', sst.type='mur', spatLim = sp.lim, save.dir = sst.dir)
  
  # YOU NEED SOME REPRESENTATION OF ENVIRONMENTAL DEPTH-TEMPERATURE
  # HYCOM DATA
  hycom.dir <- './tmp/hycom/'
  if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
  pdt.udates <- unique(as.Date(pdt$Date))
  get.env(pdt.udates, filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
  
  # BATHYMETRY
  bathy.dir <- './tmp/bathy/'
  if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
  bathy <- get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir)
  bathy <- raster(paste(bathy.dir, 'bathy.nc', sep=''))
  #library(raster); plot(bathy)
  # OR READ IT FROM NETCDF
  #bathy.nc <- RNetCDF::open.nc(paste(bathy.dir, 'bathy.nc', sep=''))

  #------------
  # CALCULATE LIKELIHOODS
  #------------
  # .par functions are the same calculations as those lacking .par, except they have been parallelized to leverage multiple CPUs
  locs.grid <- setup.locs.grid(sp.lim)
  
  # LIGHT-BASED LIKELIHOODS
  #L.light <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25) # if trying to use raw light levels, not currently recommended (v0.2)
  L.light <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
  #library(fields);library(raster)
  #plot(L.light[[12]]); world(add=T)
  
  # SST LIKELIHOODS
  L.sst <- calc.sst(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
  #L.sst <- calc.sst.par(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
  saveRDS(L.sst, file=paste0('L.sst', format(Sys.Date(), '%Y%m%d'), '.RDS')) # good idea to save after these larger calculations in case the next one causes problems
  gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
  
  # PDT LIKELIHOODS
  # OCEAN HEAT CONTENT (INTEGRATED PDTS)
  L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F, bathy = F)
  saveRDS(L.ohc, file=paste0('L.ohc', format(Sys.Date(), '%Y%m%d'), '.RDS')) # good idea to save after these larger calculations in case the next one causes problems
  gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

  # HYCOM PROFILE BASED LIKELIHOODS
  #pdt$Date <- as.POSIXct(pdt$Date, tz='UTC')
  L.hycom <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 3, dateVec = dateVec, use.se = F)
  saveRDS(L.hycom, file=paste0('L.hycom', format(Sys.Date(), '%Y%m%d'), '.RDS')) # good idea to save after these larger calculations in case the next one causes problems
  ## took ~32 mins on a tiny machine on AWS
  
  # BATHY PROFILE BASED LIKELIHOODS
  
  L.bathy <- calc.bathy(mmd, bathy.grid = bathy, dateVec, focalDim = 3, sens.err = 5)

  #save.image(file='~/work/Data/thresher/196385/196385_hmmoce_bathy_likelihoods_20200708.rda') # good idea to save after these larger calculations in case the next one causes problems
  gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

  L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.hycom, L.bathy = L.bathy)
  save(L.rasters, file=paste0(meta$instrument_name[i],'_L.rasters_', format(Sys.Date(), '%Y%m%d'), '.rda'))
  
  ## resample rasters
  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  #L.rasters$L.ohc <- raster::shift(L.rasters$L.ohc, dx = 360)
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  saveRDS(L.res, file=paste0('L.res_', format(Sys.Date(), '%Y%m%d'), '.RDS'))
}

  #------------
  # PREPARE TO RUN THE MODEL
  #------------
  run_idx <- list(c(1,2,3,5),
                  c(1,2,4,5),
                  c(1,2,5),
                  c(1,3,5),
                  c(1,4,5))
  
  ## each run idx combination
  for (tt in 1:length(run_idx)){
    
    #iniloc$lon <- make360(iniloc$lon)
    t1 <- Sys.time()
    L <- HMMoce:::make.L(L.res$L.rasters[run_idx[[tt]]], iniloc, dateVec) 
    t2 <- Sys.time()
    print(t2-t1)
    print(tt)
    
    saveRDS(L, file=paste0('L_tt', tt, '_', format(Sys.Date(), '%Y%m%d'), '.rda'))
    gc()
    
    
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
  
  
  L.4 <- L.sst * 0 # just fill the missing likelihoods (that we're currently skipping) with zeros
  for (i in 1:21) L.bathy[[i]] <- L.bathy[[i]] * 0
  #L.rasters <- mget(ls(pattern = 'L\\.')) # use with caution as all workspace items containing 'L.' will be listed. We only want the likelihood outputs calculated above
  L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.4 = L.4, L.hycom = L.hycom, L.bathy = L.bathy)
  #resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  resamp.idx <- 2
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  
  # Figure out appropriate L combinations
  # use this if you have a vector (likVec) indicating which likelihoods you are calculating
  # for example, likVec <- c(1,2,5) for light, sst, and hycom likelihoods
  likVec <- c(1,2,3,5,6) 
  if (length(likVec) > 2){
    L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
  } else{
    L.idx <- utils::combn(likVec, 2, simplify=F)
  }
  
  # which of L.idx combinations do you want to run?
  run.idx <- c(1:10)
  
  # vector of appropriate bounding in filter. see ?hmm.filter for more info
  bndVec <- c(NA, 1, 2)
  
  # vector of appropriate migr kernel speed. see ?makePar for more info.
  parVec <- c(1, 2)
  
  # GOOD IDEA TO CLEAN THINGS UP AND SAVE
  #rm(list=c('L.light','L.sst','L.ohc','L.4','L.hycom', 'woa.quarter'))
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
  ptt <- 196385
  setwd(paste('~/work/Results/thresher/',sep=''))
  
  ## 196385_idx4_bndNA_par1-HMMoce_res.rda
  tt = 18; bnd = NA; i = 1 ## tt = 18 is for likelihoods 2, 3, 6
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
        res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.hycom)[1])
        #setwd()
        save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        #save.image(file=paste(ptt, '-HMMoce.RData', sep=''))
        #source('~/HMMoce/R/hmm.diagnose.r') # not yet functional
        #hmm.diagnose(res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt, plot=T)
        
        write.table(outVec, file=paste(as.character(ptt), '_HMMoce_results_outVec.csv', sep=''), sep=',', append=T, row.names = F, col.names = F)
        
      } # parVec loop
    } # bndVec loop
  } # L.idx loop
  
  
}

#} # loop across indiv ptts
#parallel::stopCluster(cl)
#closeAllConnections()
load('~/work/Results/thresher/196385_idx4_bndNA_par1-HMMoce_res.rda')
load('~/work/Results/thresher/196385_idx7_bndNA_par2-HMMoce_res.rda')
load('~/work/Results/thresher/196385_idx18_bndNA_par1-HMMoce_res.rda')

# run getCtr
bnds <- getCtr(res$s, res$tr, res$g, threshold = 10, makePlot=F)

# compile df of position, date, bnds (xtracto)
#data.frame(lapply(bnds, FUN=function(x) c(x$yDist, x$xDist)))
df <- cbind(res$tr, t(data.frame(lapply(bnds, FUN=function(x) c(x$yDist, x$xDist)))))
names(df)[5:6] <- c('ydist','xdist')
df$ydist[which(is.na(df$ydist))] <- mean(df$ydist, na.rm=T)
df$xdist[which(is.na(df$xdist))] <- mean(df$xdist, na.rm=T)
df$ptt <- 196385
write.table(df, file='~/work/Results/thresher/196385_idx18_bndNA_par1-HMMoce_res.csv', sep=',', col.names = T, row.names = F)

#ci <- list()
ci <- CI2shp(df)
plot(df$lon, df$lat, xlim=c(37.5,39.25), ylim=c(21.3,22.8))
plot(ci, add=T, col=alpha('grey', 0.5), border=NA)
lines(df$lon, df$lat)
world(add=T)

png('196385-bathy2.png', width=12, height=10, units='in', res=300)
plot(L.bathy[[22]])
world(add=T, fill=T)
contour(L.res$g$lon[1,], L.res$g$lat[,1], L[20,,], add=T)
points(locs$Longitude, locs$Latitude)
points(iniloc$lon[1], iniloc$lat[1], pch=16, col='green')
points(iniloc$lon[2], iniloc$lat[2], pch=16, col='red')
#points(38.575, 22.125)
dev.off()

