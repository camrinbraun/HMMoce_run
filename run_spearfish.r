#========================
## HMMoce run w/example data
#========================
# might be a good idea to install latest version of HMMoce
# install.packages('HMMoce')
install_github("camrinbraun/HMMoce", ref='dev', dependencies=F)
#library(HMMoce)

#------------
# LOAD THE TAG DATA
#------------
# setwd()

# SET INITIAL LOCATIONS (TAG AND POP-UP)
meta <- read.table('~/ebs/Data/spearfish/igfa_spearfish_metadata_2018mar9.csv', sep=',', header=T)
meta <- meta[-1,]
meta$deploy_date <- as.Date(meta$deploy_date, format='%m/%d/%y')
meta$popdate <- as.Date(meta$popdate, format='%m/%d/%y')
meta <- meta[which(!is.na(meta$popdate)),]
meta$ptt <- meta$tagnumber
meta$eventid <- substr(meta$eventid,1,7)

for (zz in 1:nrow(meta)){
ptt <- meta$ptt[zz]
  iniloc <- data.frame(matrix(c(lubridate::day(meta$deploy_date[zz]), lubridate::month(meta$deploy_date[zz]), lubridate::year(meta$deploy_date[zz]), meta$deploylat[zz], meta$deploylon[zz], 
                              lubridate::day(meta$popdate[zz]), lubridate::month(meta$popdate[zz]), lubridate::year(meta$popdate[zz]), meta$poplat[zz], meta$poplon[zz]), nrow = 2, ncol = 5, byrow = T))

names(iniloc) <- list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.POSIXct(seq(tag, pop, by = '6 hours')) 

# READ IN DATA AS OUTPUT FROM WC PORTAL
# SST DATA
setwd(paste('~/ebs/Data/spearfish/tag_data/', meta$eventid[zz], '_', meta$ptt[zz],'/',sep=''))
sstFile <- paste(meta$eventid[zz],'_',meta$ptt[zz], '-SST.csv',sep='')
tag.sst <- read.wc(filename=sstFile, type = 'sst', tag=tag, pop=pop, dateFormat = '%m/%d/%Y %H:%M:%S', verbose=T) 
sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

# DEPTH-TEMPERATURE PROFILE DATA
pdtFile <- paste(meta$eventid[zz],'_',meta$ptt[zz], '-PDTs.csv',sep='')
pdt <- read.wc(pdtFile, type = 'pdt', tag=tag, pop=pop, verbose=T) 
pdt.udates <- pdt$udates; pdt <- pdt$data

seriesFile <- paste(meta$eventid[zz],'_',meta$ptt[zz], '-Series.csv',sep='')
series <- read.table(seriesFile, sep=',', header=T, blank.lines.skip = F)




# RAW LIGHT DATA
#lightFile <- system.file("extdata", "141259-LightLoc.csv", package = "HMMoce")
#light <- read.wc(ptt, lightFile, type = 'light', tag=tag, pop=pop); 
#light.udates <- light$udates; light <- light$data

# LIGHT BASED POSITIONS FROM GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
locsFile <- paste(meta$eventid[zz],'_',meta$ptt[zz], '-Locations.csv',sep='')
locs <- read.table(locsFile, sep = ',', header = T, blank.lines.skip = F)
locDates <- as.Date(as.POSIXct(locs$Date, format=findDateFormat(locs$Date)))

# SET SPATIAL LIMITS
# these are the lat/lon bounds of your study area (e.g. where you think the animal went)
sp.lim <- list(lonmin = -165,
               lonmax = -148,
               latmin = 12,
               latmax = 25)

#------------ 
##  GET ENVIRONMENTAL DATA 
#------------ 
# env data downloads can be
#large, depending on application for 180 days of data spanning the NW Atlantic
#(the example application), the downloads will take ~10mins on Amazon EC2.
#Personal computers will likely be slower.

# DOWNLOAD SST DATA
#sst.dir <- paste(tempdir(), '/sst/', sep='')
sst.dir <- '~/ebs/EnvData/sst/spearfish/'
#dir.create(sst.dir, recursive = TRUE)
get.env(sst.udates, filename='mursst', type = 'sst', sst.type='mur', spatLim = sp.lim, save.dir = sst.dir)

# YOU NEED SOME REPRESENTATION OF ENVIRONMENTAL DEPTH-TEMPERATURE
# HYCOM DATA
#hycom.dir <- paste(tempdir(), '/hycom/', sep='')
hycom.dir <- '~/ebs/EnvData/hycom/spearfish/'
#dir.create(hycom.dir, recursive = TRUE)
get.env(pdt.udates, filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

# OR WORLD OCEAN ATLAS DATA
#woa.dir <- paste(tempdir(), '/woa/', sep='')
#dir.create(woa.dir, recursive = TRUE)
#get.env(type = 'woa', resol = 'quarter', save.dir = woa.dir)
# THEN LOAD AND CHECK THE DOWNLOADED RDA FILE FOR WOA
#load(paste(woa.dir,'woa.quarter.rda',sep=''))
#str(woa.quarter)
#List of 4
#$ watertemp: num [1:44, 1:46, 1:57, 1:12] 26.5 26.5 26.4 26.3 26.2 ...
#$ lon      : num [1:44(1d)] -95.5 -94.5 -93.5 -92.5 -91.5 -90.5 -89.5 -88.5 -87.5 -86.5 ...
#$ lat      : num [1:46(1d)] 9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 ...
#$ depth    : num [1:57(1d)] 0 5 10 15 20 25 30 35 40 45 ...

# BATHYMETRY
#bathy.dir <- paste(tempdir(), '/bathy/', sep='')
bathy.dir <- '~/ebs/EnvData/bathy/spearfish/'
#dir.create(bathy.dir, recursive = TRUE)
bathy <- get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir)
#library(raster); plot(bathy)
# OR READ IT FROM NETCDF
#bathy.nc <- RNetCDF::open.nc(paste(bathy.dir, 'bathy.nc', sep=''))

#------------
# CALCULATE LIKELIHOODS
#------------
# .par functions are the same calculations as those lacking .par, except they have been parallelized to leverage multiple CPUs
locs.grid <- setup.locs.grid(sp.lim)

# vector indicating which likelihoods to run (e.g. 1=light, 2=sst, 5=hycom)
# can be combined with if() statements around calc functions: if (any(likVec == 5) & !exists('L.5')){calc.hycom(...)}
likVec <- c(1,2,5) 

# LIGHT-BASED LIKELIHOODS
#L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25) # if trying to use raw light levels, not currently recommended (v0.2)
L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
#library(fields);library(raster)
#plot(L.1[[12]]); world(add=T)

# SST LIKELIHOODS
L.2 <- calc.sst(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
#L.2 <- calc.sst.par(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
# save.image() # good idea to save after these larger calculations in case the next one causes problems
gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

# PDT LIKELIHOODS
# OCEAN HEAT CONTENT (INTEGRATED PDTS)
#L.3 <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
# save.image() # good idea to save after these larger calculations in case the next one causes problems
gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

# WORLD OCEAN ATLAS-BASED LIKELIHOODS
#L.4 <- calc.woa.par(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
# save.image() # good idea to save after these larger calculations in case the next one causes problems
gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

# HYCOM PROFILE BASED LIKELIHOODS
L.5 <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
# save.image() # good idea to save after these larger calculations in case the next one causes problems
gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
#save.image('~/ebs/example.rda')

#------------
# PREPARE TO RUN THE MODEL
#------------
L.rasters <- mget(ls(pattern = 'L\\.')) # use with caution as all workspace items containing 'L.' will be listed. We only want the likelihood outputs calculated above
resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
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
run.idx <- c(1,2,4)

# vector of appropriate bounding in filter. see ?hmm.filter for more info
bndVec <- c(NA, 5, 10)

# vector of appropriate migr kernel speed. see ?makePar for more info.
parVec <- c(2, 4)

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

for (tt in run.idx){
  for (bnd in bndVec){
    for (i in parVec){
      
      ptt=141259
      runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i,sep='')
      
      # COMBINE LIKELIHOOD MATRICES
      # L.idx combination indicates likelihood surfaces to consider
      L <- make.L(L1 = L.res[[1]][L.idx[[tt]]],
                  L.mle.res = L.res$L.mle.res, dateVec = dateVec,
                  locs.grid = locs.grid, iniloc = iniloc, bathy = bathy,
                  pdt = pdt)
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
      
      write.table(outVec, file='HMMoce_results_outVec.csv', sep=',', append=T)
      
    } # parVec loop
  } # bndVec loop
} # L.idx loop

} # loop across indiv ptts
#parallel::stopCluster(cl)
#closeAllConnections()