library(HMMoce)
#setwd('~/HMMoce/'); devtools::load_all()
dataDir <- '~/ebs/Data/BaskingSharks/batch/'
envDir <- '~/ebs/EnvData/'
sst.dir <- '~/ebs/EnvData/sst/BaskingSharks/'
hycom.dir <- '~/ebs/EnvData/hycom3/BaskingSharks/'
statusVec <- NA
bucketDir <- 'braun-data/Data/BaskingSharks/batch'

load('~/ebs/Data/BaskingSharks/batch/bask_lims.rda')
str(bask.lims)
# i think 52556, 52562, 100975 are "big" bounds. may be a few more.
# add this spec to meta

# aws credentials for reading s3
Sys.setenv("AWS_ACCESS_KEY_ID" = aws.signature::locate_credentials()[[1]],
           "AWS_SECRET_ACCESS_KEY" = aws.signature::locate_credentials()[[2]],
           "AWS_DEFAULT_REGION" = aws.signature::locate_credentials()[[4]])
library(aws.s3)
#bucketlist()

# which of L.idx combinations do you want to run?
run.idx <- c(1,2,4,7,11,13)

# vector of appropriate bounding in filter
bndVec <- c(NA, 5, 10)

# vector of appropriate migr kernel speed
parVec <- c(2, 4)

setwd(dataDir)
#aws.s3::save_object('bask_metadata.csv', file='bask_metadata.csv', bucket=bucketDir)
meta <- read.table(paste(dataDir, 'bask_metadata.csv',sep=''), sep=',', header=T)
likVec=c(1,2,3,4,5)

#for (ii in 1:nrow(meta)){ #nextAnimal
ii = 37
ptt <- meta$PTT[ii] #nextAnimal

# set an area of interest for a particular individual in the resample.grid function using:
#bnds <- raster::extent(-70, -60, 15, 33)
setwd(paste(dataDir, '/', ptt, sep=''))

# check for status. which checkpoint to start at?
tryCatch({
  err <- try(
    #aws.s3::save_object('check3.rda', file='check3.rda', bucket=paste(bucketDir, '/', ptt, sep='')),
    file.exists('check3.rda'),
    silent = T)
}, error=function(e){print(paste('ERROR: Data does not exist for this checkpoint.', sep = ''))})

if (err == FALSE){
  tryCatch({
    err <- try(
      #aws.s3::save_object('check2.rda', file='check2.rda', bucket=paste(bucketDir, '/', ptt, sep='')),
      file.exists('check2.rda'),
      silent = T)
  }, error=function(e){print(paste('ERROR: Data does not exist for this checkpoint.', sep = ''))})
  
  if (err == FALSE){
    tryCatch({
      err <- try(
        #aws.s3::save_object('check1.rda', file='check1.rda', bucket=paste(bucketDir, '/', ptt, sep='')),
        file.exists('check1.rda'),
        silent = T)
    }, error=function(e){print(paste('ERROR: Data does not exist for this checkpoint.', sep = ''))})
    
    if (err == FALSE){
      enterAt <- 1
    } else{
      enterAt <- 2
    }
    
  } else{
    enterAt <- 3
  }
  
} else{
  # skip this one and go on to next animal
  stop(paste('PTT ', ptt, ' is complete at check3. Go to ii=', ii+1, '.', sep=''))
}

statusVec <- c(paste('Entered at ', enterAt, sep='')); print(statusVec)

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
  setwd(myDir)
  #load(paste(myDir, ptt,'_likelihoods2.RData', sep=''))
  
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
  locs <- read.table(paste(myDir, ptt, '-Locations-GPE2.csv', sep = ''), sep = ',', header = T, blank.lines.skip = F)
  locDates <- as.Date(as.POSIXct(locs$Date, format=HMMoce:::findDateFormat(locs$Date)))
  
  # SET SPATIAL LIMITS
  if(meta$bnds[ii] == 'small'){
    sp.lim <- list(lonmin = bask.lims$small.bnds[1],
                   lonmax = bask.lims$small.bnds[2],
                   latmin = bask.lims$small.bnds[3],
                   latmax = bask.lims$small.bnds[4])
    #dir.create('~/ebs/EnvData/bathy/BaskingSharks/')
    #setwd('~/ebs/EnvData/bathy/BaskingSharks/')
    #aws.s3::save_object('bask_bathy_small.grd', file='bask_bathy_small.grd', bucket='gaube-data/braun/EnvData/bathy/BaskingSharks')
    #aws.s3::save_object('bask_bathy_small.gri', file='bask_bathy_small.gri', bucket='gaube-data/braun/EnvData/bathy/BaskingSharks')
    bathy <- raster::raster(paste(envDir,'bathy/BaskingSharks/bask_bathy_small.grd',sep=''))
    is.small <- TRUE
    
  } else if(meta$bnds[ii] == 'big'){
    sp.lim <- list(lonmin = bask.lims$big.bnds[1],
                   lonmax = bask.lims$big.bnds[2],
                   latmin = bask.lims$big.bnds[3],
                   latmax = bask.lims$big.bnds[4])
    bathy <- raster::raster(paste(envDir,'bathy/BaskingSharks/bask_bathy_big.grd',sep=''))
    is.small <- FALSE
  }
  
  locs.grid <- setup.locs.grid(sp.lim)
  
  # save workspace image to s3 as checkpoint
  base::save.image(file='check1.rda')
  aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check1.rda')
  
  #================
  ## END STEP 1
  #================
  
  enterAt <- 2
}

if (enterAt == 2){
  
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
    if(is.small){
      fname <- 'bask_small'
      sst.dir.small <- paste(sst.dir, 'small/', sep='')
      L.2 <- calc.sst.par(tag.sst, filename=fname, sst.dir = sst.dir.small, dateVec = dateVec, sens.err = 1)
    } else{
      fname <- 'bask_big'
      sst.dir.big <- paste(sst.dir, 'big/', sep='')
      L.2 <- calc.sst.par(tag.sst, filename=fname, sst.dir = sst.dir.big, dateVec = dateVec, sens.err = 1)
    }
    
    #raster::cellStats(L.2, 'max')
    #aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check1.rda')
    setwd(myDir); base::save.image('check1.rda')
  }
  
  if (any(likVec == 3) & !exists('L.3')){
    if(is.small){
      fname <- 'bask'
      hycom.dir <- paste(hycom.dir, 'small/', sep='')
      if(length(pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 6, 15)))]) > 0) stop('Not all hycom data is available!')
    } else{
      fname <- 'bask_big'
      hycom.dir <- paste(hycom.dir, 'big/', sep='')
      if(length(pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 10, 19)))]) > 0) stop('Not all hycom data is available!')
    }
    
    L.3 <- calc.ohc.par(pdt, filename=fname, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F, ncores=parallel::detectCores()-2)
    # checkpoint each big L calculation step
    if (exists('L.3')){
      ohc.se <- F
      #aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check1.rda')
      setwd(myDir); base::save.image('check1.rda')
    } else{
      L.3 <- calc.ohc.par(pdt, filename=fname, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = T)
      #aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check1.rda')
      setwd(myDir); base::save.image('check1.rda')
      if (!exists('L.3')){
        warning('Error: calc.ohc function failing for both standard error calculations.')
        statusVec <- c(statusVec, 'calc.ohc function failed for both standard error calculations')
      } else{
        ohc.se <- T
      }
    }
  }
  
  hycom.dir <- '~/ebs/EnvData/hycom3/BaskingSharks/'
  
  if (any(likVec == 4) & !exists('L.4')){
    load('~/ebs/EnvData/woa/woa.quarter.rda')
    L.4 <- calc.woa.par(pdt, ptt=ptt, woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = T, ncores=parallel::detectCores()-2)
    # checkpoint each big L calculation step
    if (exists('L.4')){
      woa.se <- T
      #aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check1.rda')
      setwd(myDir); base::save.image('check1.rda')
    } else{
      warning('Error: calc.woa function failed.')
      statusVec <- c(statusVec, 'calc.woa function failed')
      L.4 <- L.1 * 0
      
    }
  } else{
    L.4 <- L.1 * 0
  }
  
  if (any(likVec == 5) & !exists('L.5')){
    if(is.small){
      fname <- 'bask'
      hycom.dir <- paste(hycom.dir, 'small/', sep='')
      if(length(pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 6, 15)))]) > 0) stop('Not all hycom data is available!')
    } else{
      fname <- 'bask_big'
      hycom.dir <- paste(hycom.dir, 'big/', sep='')
      if(length(pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 10, 19)))]) > 0) stop('Not all hycom data is available!')
    }
    
    L.5 <- calc.hycom.par(pdt, filename=fname, hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
    if (exists('L.5')){
      hyc.se <- T
      #aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check1.rda')
      setwd(myDir); base::save.image('check1.rda')
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
  
  # eventually can be used to fill any empty L in L.res
  #rvec <- list()
  #for (iii in 1:5){
  #  rvec[[iii]] <- exists(paste('L.',iii,sep=''))
  #}
  
  L.rasters <- mget(ls(pattern = 'L\\.'))
  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  
  # Figure out appropriate L combinations
  if (length(likVec) > 2){
    L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
  } else{
    L.idx <- utils::combn(likVec, 2, simplify=F)
  }
  
  rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))
  
  # save workspace image to s3 as checkpoint
  setwd(myDir); base::save.image('check2.rda')
  #aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check2.rda')
  if(file.exists('check2.rda')) file.remove('check1.rda')
  
  #================
  ## END STEP 2
  #================
  enterAt <- 3
} 

if (enterAt == 3){
  
  setwd(myDir); load('check2.rda')
  
  #optim function here
  
  # BEGIN PARALLEL STUFF  
  print('Processing in parallel... ')
  ncores <- parallel::detectCores()
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  
  ans = foreach::foreach(tt = run.idx) %dopar%{
    #setwd('~/HMMoce'); devtools::load_all()
    setwd(myDir)
    library(HMMoce)
    #for (tt in run.idx){
    for (bnd in bndVec){
      for (i in parVec){
        
        runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i,sep='')
        
        #----------------------------------------------------------------------------------#
        # COMBINE LIKELIHOOD MATRICES
        #----------------------------------------------------------------------------------#
        L <- HMMoce::make.L(L1 = L.res[[1]][L.idx[[tt]]],
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
        par0 <- HMMoce::makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
        #K1 <- par0$K1; K2 <- par0$K2; 
        P.final <- par0$P.final
        
        # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
        par0 <- HMMoce::makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
        K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
        
        # RUN THE FILTER STEP
        if(!is.na(bnd)){
          f <- HMMoce::hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
          maskL.logical <- TRUE
        } else{
          f <- HMMoce::hmm.filter(g, L, K1, K2, P.final, maskL=F)
          maskL.logical <- FALSE
        }
        nllf <- -sum(log(f$psi[f$psi>0]))
        
        # RUN THE SMOOTHING STEP
        s <- HMMoce::hmm.smoother(f, K1, K2, L, P.final)
        
        # GET THE MOST PROBABLE TRACK
        tr <- HMMoce::calc.track(s, g, dateVec)
        setwd(myDir)
        
        # WRITE OUT RESULTS
        outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
                           Lidx = paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                           spLims = sp.lim[1:4], resol = raster::res(L.rasters[[resamp.idx]]),
                           maskL = maskL.logical, NLL = nllf, name = runName), ncol=15)
        #write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
        #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
        res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
        setwd(myDir); save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        #save.image(file=paste(ptt, '-HMMoce.RData', sep=''))
        aws.s3::s3save(res, bucket=paste(bucketDir, '/', ptt, sep=''), object=paste(runName, '-HMMoce_res.rda', sep=''))
        source('~/HMMoce/R/hmm.diagnose.r')
        hmm.diagnose(res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt, plot=T)
        
        outVec <- outVec
        
      } # parVec loop
    } # bndVec loop
  } # L.idx loop
  
  
  parallel::stopCluster(cl)
  
}

# save workspace image to s3 as checkpoint
#setwd(myDir); base::save.image('check3.rda')
#aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check3.rda')


