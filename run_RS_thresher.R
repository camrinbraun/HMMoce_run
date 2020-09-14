#========================
## HMMoce run w/example data
#========================
# might be a good idea to install latest version of HMMoce
# install.packages('HMMoce')
#devtools::install_github("camrinbraun/HMMoce", ref='dev', dependencies=F)
#library(HMMoce); 
devtools::load_all('~/work/RCode/HMMoce')
library(raster)

# get necessary funs
#install.packages(c('readxl','readr')) #lubridate, plyr, MASS
library(readr); library(readxl); library(lubridate)
#devtools::install_github('camrinbraun/HMMoce', dependencies = F)
#source('https://raw.githubusercontent.com/galuardi/analyzepsat/v4.0/R/MWTextract.r')

## load all tag metadata master sheet
meta <- read.table('~/work/RCode/nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, blank.lines.skip = F, skip=0)


#------------
# LOAD THE TAG DATA
#------------
# setwd()

# SET INITIAL LOCATIONS (TAG AND POP-UP)
## we currently use a unique instrument name for each deployment (taxonserialnumber_deploymentyear_tagPTT)
idx <- which(meta$instrument_name == '159922_2020_196385')
tag <- as.POSIXct(meta$time_coverage_start[idx], format='%Y-%m-%d %H:%M:%S', tz='UTC')
pop <- as.POSIXct(meta$time_coverage_end[idx], format='%Y-%m-%d %H:%M:%S', tz='UTC')

#iniloc <- data.frame(matrix(c(lubridate::day(tag), lubridate::month(tag), lubridate::year(tag),
#                              meta$geospatial_lat_start[idx], meta$geospatial_lon_start[idx], 
#                              lubridate::day(pop), lubridate::month(pop), lubridate::year(pop),
#                              meta$geospatial_lat_end[idx], meta$geospatial_lon_end[idx]),
#                            nrow = 2, ncol = 5, byrow = T))
iniloc <- data.frame(matrix(c(lubridate::day(tag), lubridate::month(tag), lubridate::year(tag),
                              meta$geospatial_lat_start[idx], meta$geospatial_lon_start[idx], 
                              lubridate::day(pop), lubridate::month(pop), lubridate::year(pop),
                              22.125, 38.575),
                            nrow = 2, ncol = 5, byrow = T))
names(iniloc) <- list('day','month','year','lat','lon')
#tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
#pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA AS OUTPUT FROM WC PORTAL
# SST DATA
sstFile <- '~/work/Data/thresher/196385/196385-SST.csv'
tag.sst <- read.wc(sstFile, type = 'sst', tag=tag, pop=pop, verbose=T, dateFormat = '%H:%M:%S %d-%b-%Y') 
sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

# DEPTH-TEMPERATURE PROFILE DATA
pdtFile <- '~/work/Data/thresher/196385/196385-PDTs.csv'
pdt <- read.wc(pdtFile, type = 'pdt', tag=tag, pop=pop, verbose=T) 
pdt.udates <- pdt$udates; pdt <- pdt$data

# RAW LIGHT DATA
#lightFile <- system.file("extdata", "141259-LightLoc.csv", package = "HMMoce")
#light <- read.wc(ptt, lightFile, type = 'light', tag=tag, pop=pop); 
#light.udates <- light$udates; light <- light$data

# LIGHT BASED POSITIONS FROM GPE2 (INSTEAD OF RAW LIGHTLOCS FROM PREVIOUS)
locsFile <- '~/work/Data/thresher/196385/196385-Locations-GPE2.csv'
locs <- read.table(locsFile, sep = ',', header = T, blank.lines.skip = F)
locDates <- as.Date(as.POSIXct(locs$Date, format=findDateFormat(locs$Date)))

# SET SPATIAL LIMITS
# these are the lat/lon bounds of your study area (e.g. where you think the animal went)
sp.lim <- list(lonmin = 37.5,
               lonmax = 39,
               latmin = 21.5,
               latmax = 22.5)

#------------ 
##  GET ENVIRONMENTAL DATA 
#------------ 
# env data downloads can be
#large, depending on application for 180 days of data spanning the NW Atlantic
#(the example application), the downloads will take ~10mins on Amazon EC2.
#Personal computers will likely be slower.

# DOWNLOAD SST DATA
#sst.dir <- paste(tempdir(), '/sst/', sep='')
sst.dir <- '~/work/EnvData/sst/thresher/'
dir.create(sst.dir, recursive = TRUE)
sst.udates <- unique(as.Date(tag.sst$Date))
get.env(sst.udates, filename='mursst', type = 'sst', sst.type='mur', spatLim = sp.lim, save.dir = sst.dir)

# YOU NEED SOME REPRESENTATION OF ENVIRONMENTAL DEPTH-TEMPERATURE
# HYCOM DATA
#hycom.dir <- paste(tempdir(), '/hycom/', sep='')
hycom.dir <- '~/work/EnvData/hycom/thresher/'
dir.create(hycom.dir, recursive = TRUE)
pdt.udates <- unique(as.Date(pdt$Date))
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
bathy.dir <- '~/work/EnvData/bathy/thresher/'
#dir.create(bathy.dir, recursive = TRUE)
#bathy <- get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir)
#bathy <- raster(paste(bathy.dir, 'bathy.nc', sep=''))
#library(raster); plot(bathy)
# OR READ IT FROM NETCDF
#bathy.nc <- RNetCDF::open.nc(paste(bathy.dir, 'bathy.nc', sep=''))
bathy <- raster('~/work/EnvData/bathy/red_sea_swords.grd')

#------------
# CALCULATE LIKELIHOODS
#------------
# .par functions are the same calculations as those lacking .par, except they have been parallelized to leverage multiple CPUs
locs.grid <- setup.locs.grid(sp.lim)

# vector indicating which likelihoods to run (e.g. 1=light, 2=sst, 5=hycom)
# can be combined with if() statements around calc functions: if (any(likVec == 5) & !exists('L.5')){calc.hycom(...)}
likVec <- c(1,2,3,5,6) 

# LIGHT-BASED LIKELIHOODS
#L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25) # if trying to use raw light levels, not currently recommended (v0.2)
L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
#library(fields);library(raster)
#plot(L.1[[12]]); world(add=T)

# SST LIKELIHOODS
L.2 <- calc.sst(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
#L.2 <- calc.sst.par(tag.sst, filename='mursst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
save.image(file='~/ebs/Data/AndersonJ/161211/161211_l2.rda') # good idea to save after these larger calculations in case the next one causes problems
gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

# PDT LIKELIHOODS
# OCEAN HEAT CONTENT (INTEGRATED PDTS)
L.3 <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F, bathy = F)
# save.image() # good idea to save after these larger calculations in case the next one causes problems
gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

# WORLD OCEAN ATLAS-BASED LIKELIHOODS
#L.4 <- calc.woa.par(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
# save.image() # good idea to save after these larger calculations in case the next one causes problems
gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running

# HYCOM PROFILE BASED LIKELIHOODS
#pdt$Date <- as.POSIXct(pdt$Date, tz='UTC')
L.5 <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 3, dateVec = dateVec, use.se = F)
## took ~32 mins on a tiny machine on AWS

mmd <- read.table('~/work/Data/thresher/196385/196385-MinMaxDepth.csv', sep=',', header=T)
mmd$Date <- as.Date(as.POSIXct(mmd$Date, format='%H:%M:%S %d-%b-%Y', tz='UTC'))
mmd <- mmd[which(!is.na(mmd$MaxDepth)),]
mmd <- mmd[which(!duplicated(mmd$Date)),]
L.6 <- calc.bathy(mmd, bathy.grid = bathy, dateVec, focalDim = 3, sens.err = 5)
#calc.bathy=function (tag.pdt,bathy.grid, dateVec, focalDim = NULL, sens.err = 1){
  
save.image(file='~/work/Data/thresher/196385/196385_hmmoce_bathy_likelihoods_20200708.rda') # good idea to save after these larger calculations in case the next one causes problems
#gc(); closeAllConnections() # also good to do garbage collection and kill any straggling processes that are running
#save.image('~/ebs/Data/AndersonJ/161211_Likelihoods.rda')

#------------
# PREPARE TO RUN THE MODEL
#------------

#load(paste('~/ebs/Data/AndersonJ/161211_Likelihoods.rda',sep=''))
#known$lat <- as.numeric(as.character(known$lat))
#known$lon <- as.numeric(as.character(known$lon))

L.4 <- L.2 * 0 # just fill the missing likelihoods (that we're currently skipping) with zeros
for (i in 1:21) L.6[[i]] <- L.6[[i]] * 0
#L.rasters <- mget(ls(pattern = 'L\\.')) # use with caution as all workspace items containing 'L.' will be listed. We only want the likelihood outputs calculated above
L.rasters <- list(L.1 = L.1, L.2 = L.2, L.3 = L.3, L.4 = L.4, L.5 = L.5, L.6 = L.6)
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
plot(L.6[[22]])
world(add=T, fill=T)
contour(L.res$g$lon[1,], L.res$g$lat[,1], L[20,,], add=T)
points(locs$Longitude, locs$Latitude)
points(iniloc$lon[1], iniloc$lat[1], pch=16, col='green')
points(iniloc$lon[2], iniloc$lat[2], pch=16, col='red')
#points(38.575, 22.125)
dev.off()

