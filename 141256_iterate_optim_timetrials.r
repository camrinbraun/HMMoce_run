load('141256_likelihoods.RData')

t005 <- Sys.time()
rm(L.4); rm(L.idx); rm(L.rasters); rm(L.res)
L.rasters <- mget(ls(pattern = 'L\\.'))
#resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 4), mle.res=1)
L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 3), mle.res=1)
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 2), mle.res=1)
#L.res <- resample.grid(L.rasters, L.rasters[[1]], mle.res=1)
#L.res <- resample.grid(L.rasters, L.rasters[[3]], mle.res=1)


# Figure out appropriate L combinations
if (length(likVec) > 2){
  L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  L.idx <- utils::combn(likVec, 2, simplify=F)
}

rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))

bndVec <- c(NA, 5, 10)
parVec <- c(2, 4)
run.idx <- c(1,2,4,7,11,13)
require(foreach)
print('Processing in parallel... ')
ncores <- ceiling(parallel::detectCores() * .25)
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl, cores = ncores)

#optim function here
#tt <- 13
ans = foreach::foreach(tt = run.idx) %dopar%{
  
  library(HMMoce)
   
  for (bnd in bndVec){
    for (i in parVec){
      
      runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i, '_res', round(raster::res(L.res[[1]][[1]])[1],2), sep='')
      
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
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
      par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
      #K1 <- par0$K1; K2 <- par0$K2; 
      P.final <- par0$P.final
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
      par0 <- makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
      K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
      
      # RUN THE FILTER STEP
      if(!is.na(bnd)){
        f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
        maskL.logical <- TRUE
      } else{
        f <- hmm.filter(g, L, K1, K2, P.final, maskL=F)
        maskL.logical <- FALSE
      }
      nllf <- -sum(log(f$psi[f$psi>0]))
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K1, K2, L, P.final)
      
      # GET THE MOST PROBABLE TRACK
      tr <- calc.track(s, g, dateVec, iniloc)
      #setwd(myDir); plotHMM(s, tr, dateVec, ptt=runName, save.plot = T)
      
      
      # WRITE OUT RESULTS
      outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
                         Lidx = paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                         spLims = sp.lim[1:4], resol = raster::res(L.rasters[[resamp.idx]]),
                         maskL = maskL.logical, NLL = nllf, name = runName), ncol=15)
      #write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
      #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
      res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
      #setwd(myDir); 
      if(tt == 1 & i == 2){
        save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        
      }
      
    } # parVec loop
  } # bndVec loop
} # L.idx loop


parallel::stopCluster(cl)
closeAllConnections()

t006 <- Sys.time()
#===============================


load('141256_likelihoods.RData')

t007 <- Sys.time()
rm(L.4); rm(L.idx); rm(L.rasters); rm(L.res)
L.rasters <- mget(ls(pattern = 'L\\.'))
#resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 4), mle.res=1)
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 3), mle.res=1)
L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 2), mle.res=1)
#L.res <- resample.grid(L.rasters, L.rasters[[1]], mle.res=1)
#L.res <- resample.grid(L.rasters, L.rasters[[3]], mle.res=1)


# Figure out appropriate L combinations
if (length(likVec) > 2){
  L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  L.idx <- utils::combn(likVec, 2, simplify=F)
}

rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))

bndVec <- c(NA, 5, 10)
parVec <- c(2, 4)
run.idx <- c(1,2,4,7,11,13)
require(foreach)
print('Processing in parallel... ')
ncores <- ceiling(parallel::detectCores() * .25)
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl, cores = ncores)

#optim function here
#tt <- 13
ans = foreach::foreach(tt = run.idx) %dopar%{
  
  library(HMMoce)
  
  for (bnd in bndVec){
    for (i in parVec){
      
      runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i, '_res', round(raster::res(L.res[[1]][[1]])[1],2), sep='')
      
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
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
      par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
      #K1 <- par0$K1; K2 <- par0$K2; 
      P.final <- par0$P.final
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
      par0 <- makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
      K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
      
      # RUN THE FILTER STEP
      if(!is.na(bnd)){
        f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
        maskL.logical <- TRUE
      } else{
        f <- hmm.filter(g, L, K1, K2, P.final, maskL=F)
        maskL.logical <- FALSE
      }
      nllf <- -sum(log(f$psi[f$psi>0]))
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K1, K2, L, P.final)
      
      # GET THE MOST PROBABLE TRACK
      tr <- calc.track(s, g, dateVec, iniloc)
      #setwd(myDir); plotHMM(s, tr, dateVec, ptt=runName, save.plot = T)
      
      
      # WRITE OUT RESULTS
      outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
                         Lidx = paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                         spLims = sp.lim[1:4], resol = raster::res(L.rasters[[resamp.idx]]),
                         maskL = maskL.logical, NLL = nllf, name = runName), ncol=15)
      #write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
      #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
      res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
      #setwd(myDir); 
      if(tt == 1 & i == 2){
        save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        
      }      
    } # parVec loop
  } # bndVec loop
} # L.idx loop


parallel::stopCluster(cl)
closeAllConnections()



t008 <- Sys.time()

#===============================


load('141256_likelihoods.RData')

t009 <- Sys.time()
rm(L.4); rm(L.idx); rm(L.rasters); rm(L.res)
L.rasters <- mget(ls(pattern = 'L\\.'))
#resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 4), mle.res=1)
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 3), mle.res=1)
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 2), mle.res=1)
L.res <- resample.grid(L.rasters, L.rasters[[1]], mle.res=1)
#L.res <- resample.grid(L.rasters, L.rasters[[3]], mle.res=1)


# Figure out appropriate L combinations
if (length(likVec) > 2){
  L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  L.idx <- utils::combn(likVec, 2, simplify=F)
}

rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))

bndVec <- c(NA, 5, 10)
parVec <- c(2, 4)
run.idx <- c(1,2,4,7,11,13)
require(foreach)
print('Processing in parallel... ')
ncores <- ceiling(parallel::detectCores() * .25)
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl, cores = ncores)

#optim function here
#tt <- 13
ans = foreach::foreach(tt = run.idx) %dopar%{
  
  library(HMMoce)
  
  for (bnd in bndVec){
    for (i in parVec){
      
      runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i, '_res', round(raster::res(L.res[[1]][[1]])[1],2), sep='')
      
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
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
      par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
      #K1 <- par0$K1; K2 <- par0$K2; 
      P.final <- par0$P.final
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
      par0 <- makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
      K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
      
      # RUN THE FILTER STEP
      if(!is.na(bnd)){
        f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
        maskL.logical <- TRUE
      } else{
        f <- hmm.filter(g, L, K1, K2, P.final, maskL=F)
        maskL.logical <- FALSE
      }
      nllf <- -sum(log(f$psi[f$psi>0]))
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K1, K2, L, P.final)
      
      # GET THE MOST PROBABLE TRACK
      tr <- calc.track(s, g, dateVec, iniloc)
      #setwd(myDir); plotHMM(s, tr, dateVec, ptt=runName, save.plot = T)
      
      
      # WRITE OUT RESULTS
      outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
                         Lidx = paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                         spLims = sp.lim[1:4], resol = raster::res(L.rasters[[resamp.idx]]),
                         maskL = maskL.logical, NLL = nllf, name = runName), ncol=15)
      #write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
      #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
      res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
      #setwd(myDir); 
      if(tt == 1 & i == 2){
        save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        
      }      
    } # parVec loop
  } # bndVec loop
} # L.idx loop


parallel::stopCluster(cl)
closeAllConnections()


t010 <- Sys.time()

#===============================


load('141256_likelihoods.RData')

t011 <- Sys.time()
rm(L.4); rm(L.idx); rm(L.rasters); rm(L.res)
L.rasters <- mget(ls(pattern = 'L\\.'))
#resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 4), mle.res=1)
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 3), mle.res=1)
#L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 2), mle.res=1)
#L.res <- resample.grid(L.rasters, L.rasters[[1]], mle.res=1)
L.res <- resample.grid(L.rasters, L.rasters[[3]], mle.res=1)


# Figure out appropriate L combinations
if (length(likVec) > 2){
  L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  L.idx <- utils::combn(likVec, 2, simplify=F)
}

rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))

bndVec <- c(NA, 5, 10)
parVec <- c(2, 4)
run.idx <- c(1,2,4,7,11,13)
require(foreach)
print('Processing in parallel... ')
ncores <- ceiling(parallel::detectCores() * .25)
cl = parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl, cores = ncores)

#optim function here
#tt <- 13
ans = foreach::foreach(tt = run.idx) %dopar%{
  
  library(HMMoce)
  
  for (bnd in bndVec){
    for (i in parVec){
      
      runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i, '_res', round(raster::res(L.res[[1]][[1]])[1],2), sep='')
      
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
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
      par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
      #K1 <- par0$K1; K2 <- par0$K2; 
      P.final <- par0$P.final
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
      par0 <- makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
      K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
      
      # RUN THE FILTER STEP
      if(!is.na(bnd)){
        f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
        maskL.logical <- TRUE
      } else{
        f <- hmm.filter(g, L, K1, K2, P.final, maskL=F)
        maskL.logical <- FALSE
      }
      nllf <- -sum(log(f$psi[f$psi>0]))
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K1, K2, L, P.final)
      
      # GET THE MOST PROBABLE TRACK
      tr <- calc.track(s, g, dateVec, iniloc)
      #setwd(myDir); plotHMM(s, tr, dateVec, ptt=runName, save.plot = T)
      
      
      # WRITE OUT RESULTS
      outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i,
                         Lidx = paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2],
                         spLims = sp.lim[1:4], resol = raster::res(L.rasters[[resamp.idx]]),
                         maskL = maskL.logical, NLL = nllf, name = runName), ncol=15)
      #write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
      #names(outVec) <- c('ptt','bnd','migr.spd','Lidx','P1','P2','spLims','resol','maskL','nll','name')
      res <- list(outVec = outVec, s = s, g = g, tr = tr, dateVec = dateVec, iniloc = iniloc, grid = raster::res(L.res[[1]]$L.5)[1])
      #setwd(myDir); 
      if(tt == 1 & i == 2){
        save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        
      }      
    } # parVec loop
  } # bndVec loop
} # L.idx loop


parallel::stopCluster(cl)
closeAllConnections()



t012 <- Sys.time()

