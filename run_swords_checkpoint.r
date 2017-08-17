#library(HMMoce)
setwd('~/HMMoce/'); devtools::load_all()
dataDir <- '~/Data/Swordfish/batch/'
envDir <- '~/EnvData/'
sst.dir <- '~/EnvData/sst/Swordfish/'
hycom.dir <- '~/EnvData/hycom3/Swordfish/'
statusVec <- NA

# which of L.idx combinations do you want to run?
run.idx <- c(1:4, 11:16)

# vector of appropriate bounding in filter
bndVec <- c(NA, 5, 10)

# vector of appropriate migr kernel speed
parVec <- c(1, 2, 4)

meta <- read.table(paste(dataDir, 'swords_meta.csv',sep=''), sep=',', header=T)
likVec=c(1,2,3,5)

#for (ii in 1:nrow(meta)){ #nextAnimal

ptt <- meta$PTT[ii] #nextAnimal

bucketDir <- 'gaube-data/braun/Data/Swordfish/batch/'

maskL.logical <- TRUE

# check for status. which checkpoint to start at?
tryCatch({
  err <- try(
    aws.s3::save_object('check3.rda', file='check3.rda', bucket=paste(bucketDir, ptt, sep='')),
    silent = T)
}, error=function(e){print(paste('ERROR: Data does not exist for this checkpoint.', sep = ''))})

if (class(err) == 'try-error'){
  tryCatch({
    err <- try(
      aws.s3::save_object('check2.rda', file='check2.rda', bucket=paste(bucketDir, ptt, sep='')),
      silent = T)
  }, error=function(e){print(paste('ERROR: Data does not exist for this checkpoint.', sep = ''))})
  
  if (class(err) == 'try-error'){
    tryCatch({
      err <- try(
        aws.s3::save_object('check1.rda', file='check1.rda', bucket=paste(bucketDir, ptt, sep='')),
        silent = T)
    }, error=function(e){print(paste('ERROR: Data does not exist for this checkpoint.', sep = ''))})
    
    if (class(err) == 'try-error'){
      enterAt <- 1
    } else{
      enterAt <- 2
    }
    
  } else{
    enterAt <- 3
  }
  
} else{
  # skip this one and go on to next animal
  nextAnimal
}

statusVec <- c(paste('Entered at ', enterAt, sep=''))

  if (enterAt == 1){
    
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
    #bathy <- raster::raster('~/EnvData/bathy/Swordfish/sword_bathy.grd')
    aws.s3::save_object('sword_bathy.rda', file='sword_bathy.rda', bucket=paste(bucketDir, sep=''))
    load('sword_bathy.rda')
    
    # save workspace image to s3 as checkpoint
    aws.s3::s3save_image(bucket=paste(bucketDir, ptt, sep=''), object='check1.rda')
    
    #================
    ## END STEP 1
    #================
    
  } else if (enterAt == 2){
    
    load('check1.rda')
    
    #----------------------------------------------------------------------------------#
    # CALCULATE ALL LIKELIHOODS
    #----------------------------------------------------------------------------------#
    if (any(likVec == 1) & !exists('L.1')){
      #L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25)
      L.1 <- calc.gpe2(locs, locDates, iniloc = iniloc, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
      #raster::cellStats(L.1, 'max')
    }
    
    if (any(likVec == 2) & !exists('L.2')){
      L.2 <- calc.sst.par(tag.sst, filename='sword', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
      #raster::cellStats(L.2, 'max')
    }
    
    if (any(likVec == 3) & !exists('L.3')){
      if(length(pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 7, 16)))]) > 0) stop('Not all hycom data is available!')
      L.3 <- calc.ohc(pdt, filename='sword', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
      # checkpoint each big L calculation step
      if (exists('L.3')){
        ohc.se <- F
        aws.s3::s3save_image(bucket=paste(bucketDir, ptt, sep=''), object='check2.rda')
      } else{
        L.3 <- calc.ohc(pdt, filename='sword', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = T)
        aws.s3::s3save_image(bucket=paste(bucketDir, ptt, sep=''), object='check2.rda')
        if (!exists('L.3')){
          warning('Error: calc.ohc function failing for both standard error calculations.')
          statusVec <- c(statusVec, 'calc.ohc function failed for both standard error calculations')
        } else{
          ohc.se <- T
        }
      }
    }
    
    if (any(likVec == 4) & !exists('L.4')){
      L.4 <- calc.woa.par(pdt, filename='', woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = T)
      # checkpoint each big L calculation step
      if (exists('L.4')){
        woa.se <- T
        aws.s3::s3save_image(bucket=paste(bucketDir, ptt, sep=''), object='check2.rda')
      } else{
        warning('Error: calc.woa function failed.')
        statusVec <- c(statusVec, 'calc.woa function failed')
      }
    }
    
    if (any(likVec == 5) & !exists('L.5')){
      L.5 <- calc.hycom.par(pdt, filename='sword', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
      if (exists('L.5')){
        hyc.se <- T
        aws.s3::s3save_image(bucket=paste(bucketDir, ptt, sep=''), object='check2.rda')
      } else{
        warning('Error: calc.hycom function failed.')
        statusVec <- c(statusVec, 'calc.hycom.par failed')
      }
    }
    
    #----------------------------------------------------------------------------------#
    # LIST, RESAMPLE, SAVE
    #----------------------------------------------------------------------------------#
    for (iii in likVec){
      if(!exists(paste('L.',iii,sep=''))){
        statusVec <- c(statusVec, paste('Error: L.', iii, ' does not exist. Cannot combine required likelihoods.',sep=''))
        stop(paste('Error: L.', iii, ' does not exist. Cannot combine required likelihoods.',sep=''))
      } 
    }
    
    L.rasters <- mget(ls(pattern = 'L\\.'))
    resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
    L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
    
    # Figure out appropriate L combinations
    if (length(likVec) > 2){
      L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
    } else{
      L.idx <- utils::combn(likVec, 2, simplify=F)
    }
    
    # save workspace image to s3 as checkpoint
    aws.s3::s3save_image(bucket=paste(bucketDir, ptt, sep=''), object='check2.rda')
    
    #================
    ## END STEP 2
    #================
    
  } else if (enterAt == 3){
    
    load('check2.rda')
    
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
            f <- hmm.filter.ext(g, L, K1, K2, maskL=maskL.logical, P.final, minBounds = bnd)
          } else{
            f <- hmm.filter(g, L, K1, K2, P.final)
          }
          
          # RUN THE SMOOTHING STEP
          s <- hmm.smoother(f, K1, K2, L, P.final)
          
          # GET THE MOST PROBABLE TRACK
          tr <- calc.track(s, g, dateVec)
          #plotHMM(s, tr, dateVec, ptt, save.plot = F)
          
          # WRITE OUT RESULTS
          outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
                             paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                             sp.lim[1:4], raster::res(L.rasters[[resamp.idx]]),
                             L.idx[[tt]], maskL.logical), ncol=30)
          write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
          
          # checkpoint for this loop
          # fname includes paste(par, bnd, idx)
          aws.s3::s3save(a, bucket='gaube-data/braun/try_s3', object='check1.rda')
          
        } # parVec loop
      } # bndVec loop
    } # L.idx loop
    
  }
  
# save workspace image to s3 as checkpoint
aws.s3::s3save_image(bucket=paste(bucketDir, ptt, sep=''), object='check3.rda')







