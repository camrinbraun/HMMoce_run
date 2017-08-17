#library(HMMoce)
setwd('~/HMMoce/'); devtools::load_all()
dataDir <- '~/Data/Swordfish/batch/'
envDir <- '~/EnvData/'

meta <- read.table(paste(dataDir, 'swords_meta.csv',sep=''), sep=',', header=T)

likVec=c(1,2,3,5)

for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  iniloc <- data.frame(matrix(c(meta$TagDay[ii], meta$TagMonth[ii], meta$TagYear[ii], meta$TagLat[ii], meta$TagLong[ii],
                                meta$PopDay[ii], meta$PopMonth[ii], meta$PopYear[ii], meta$PopLat[ii], meta$PopLong[ii]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')

  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day'))

  # READ IN DATA FROM WC FILES
  myDir <- paste(dataDir, ptt, '/', sep='')
  load(paste(myDir, ptt,'_likelihoods2.RData', sep=''))

  # sst data
  tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop);
  sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

  # depth-temp profile data
  pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop);
  pdt.udates <- pdt$udates; pdt <- pdt$data

  # light data
  #light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop);
  #light.udates <- light$udates; light <- light$data

  # use GPE2 locs for light-based likelihoods
  locs <- read.table(paste(myDir, ptt, '-Locations.csv', sep = ''), sep = ',', header = T, blank.lines.skip = F)
  locDates <- as.Date(as.POSIXct(locs$Date, format=HMMoce:::findDateFormat(locs$Date)))

  # SET SPATIAL LIMITS
  sp.lim <- list(lonmin=-85, lonmax=-15,
                 latmin=8, latmax=53)
  locs.grid <- setup.locs.grid(sp.lim)

  # GET BATHYMETRY
  #bathy <- get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, res = c(.5))
  #raster::writeRaster(bathy, '~/EnvData/bathy/Swordfish/sword_bathy.grd')
  bathy <- raster::raster('~/EnvData/bathy/Swordfish/sword_bathy.grd')

  #----------------------------------------------------------------------------------#
  # CALCULATE ALL LIKELIHOODS
  #----------------------------------------------------------------------------------#
  if (any(likVec == 1)){
    #L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25)
    L.1 <- calc.gpe2(locs, locDates, iniloc = iniloc, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
    raster::cellStats(L.1, 'max')
  }

  if (any(likVec == 2)){
    sst.dir <- '~/EnvData/sst/Swordfish/'
    L.2 <- calc.sst.par(tag.sst, filename='sword', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
    raster::cellStats(L.2, 'max')
  }

  if (any(likVec == 3)){
    hycom.dir <- '~/EnvData/hycom3/Swordfish/'
    if(length(pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 7, 16)))]) > 0) stop('Not all hycom data is available!')
    L.3 <- calc.ohc(pdt, filename='sword', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  }

  if (any(likVec == 4)){
    L.4 <- calc.woa.par(pdt, filename='', woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = T)
  }

  if (any(likVec == 5)){
    L.5 <- calc.hycom.par(pdt, filename='sword', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
  }

  L.1 <- L.2 * 0; L.4 <- L.2 * 0

  #----------------------------------------------------------------------------------#
  # LIST, RESAMPLE, SAVE
  #----------------------------------------------------------------------------------#

  L.rasters <- mget(ls(pattern = 'L\\.'))
  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  save.image(paste(myDir, ptt, '_likelihoods3.RData', sep=''))

  # Figure out appropriate L combinations
  if (length(likVec) > 2){
    L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
  } else{
    L.idx <- utils::combn(likVec, 2, simplify=F)
  }
  run.idx <- c(1:4, 11:16)

  for (tt in run.idx){

    #----------------------------------------------------------------------------------#
    # COMBINE LIKELIHOOD MATRICES
    #----------------------------------------------------------------------------------#
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

    #bnd <- 10
    for (bnd in bndVec){
      for (i in parVec){

        # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
        par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
        #K1 <- par0$K1; K2 <- par0$K2;
        P.final <- par0$P.final

        # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
        par0 <- makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
        K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final

        # RUN THE FILTER STEP
        if(!is.na(bnd)){
          f <- hmm.filter.ext(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
        } else{
          f <- hmm.filter(g, L, K1, K2, P.final)
        }

        # RUN THE SMOOTHING STEP
        s <- hmm.smoother(f, K1, K2, L, P.final)

        # GET THE MOST PROBABLE TRACK
        tr <- calc.track(s, g, dateVec)
        #plotHMM(s, tr, dateVec, ptt, save.plot = F)

        # COMPARE HMM, GPE3, SPOT
        setwd(myDir)
        write.table(tr, file=paste(ptt, '_HMM_track.csv', sep=''), sep = ',', col.names = TRUE)
        #save.image(file=paste(ptt, '_finaltrack.RData', sep=''))

      } # parVec loop
    } # bndVec loop
  } # L.idx loop


}

#=================================
save.image(paste(myDir, ptt, '_likelihoods2.RData', sep=''))

rm(list=ls()); gc(); closeAllConnections()

#library(HMMoce)
setwd('~/HMMoce/'); devtools::load_all()
dataDir <- '~/Data/Swordfish/batch/'
envDir <- '~/EnvData/'

meta <- read.table(paste(dataDir, 'swords_meta.csv',sep=''), sep=',', header=T)

likVec=c(1,2,3,5)

for (ii in 17:nrow(meta)){
  setwd('~/HMMoce/'); devtools::load_all()
  dataDir <- '~/Data/Swordfish/batch/'
  envDir <- '~/EnvData/'

  meta <- read.table(paste(dataDir, 'swords_meta.csv',sep=''), sep=',', header=T)

  likVec=c(1,2,3,5)

  ptt <- meta$PTT[ii]
  iniloc <- data.frame(matrix(c(meta$TagDay[ii], meta$TagMonth[ii], meta$TagYear[ii], meta$TagLat[ii], meta$TagLong[ii],
                                meta$PopDay[ii], meta$PopMonth[ii], meta$PopYear[ii], meta$PopLat[ii], meta$PopLong[ii]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')

  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day'))

  # READ IN DATA FROM WC FILES
  myDir <- paste(dataDir, ptt, '/', sep='')
  load(paste(myDir, ptt,'_likelihoods.RData', sep=''))

  # depth-temp profile data
  pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop);
  pdt.udates <- pdt$udates; pdt <- pdt$data

  if (any(likVec == 3)){
    hycom.dir <- '~/EnvData/hycom3/Swordfish/'
    if(length(pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 7, 16)))]) > 0) stop('Not all hycom data is available!')
    #L.3t <- calc.ohc.par(pdt, filename='sword', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
    L.3f <- calc.ohc(pdt, filename='sword', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  }

  if (any(likVec == 5)){
    L.5t <- calc.hycom(pdt, filename='sword', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
    #L.5f <- calc.hycom.par(pdt, filename='sword', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = F)
  }
  save.image(paste(myDir, ptt, '_likelihoods2.RData', sep=''))
  res <- list(L.3 = L.3f, L.5 = L.5t)
  save(res, file=paste(myDir, ptt, '_likelihoods2_res.rda', sep=''))
  rm(L.3f); rm(L.5t); rm(pdt); rm(ptt); rm(res); gc(); closeAllConnections()
}

# 106795 failed calc.hycom.par => "newsplit: out of vertex space"

