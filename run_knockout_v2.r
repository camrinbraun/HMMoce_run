#devtools::load_all()
dataDir <- '/home/rstudio/HMMoce_run/data/'
ptt <- 141259
myDir <- paste(dataDir, ptt, '/', sep='')
setwd(myDir)
load(paste(dataDir, ptt, '/', ptt,'_finaltrack.RData', sep=''))
#load(paste(dataDir, ptt, '/', ptt,'_knock_track.RData', sep=''))
gpeNo
res$mean.gcd[5:6]
res$rmse.lat[5:6]
res$rmse.lon[5:6]

## change resample grid size???
raster::res(L.res[[1]]$L.5)

#gpeNo <- 4 #257
#gpeNo <- 5 #256
gpeNo <- 7 #254
#gpeNo <- 8 #259

source('/home/rstudio/HMMoce_run/formatTracks.r')
source('/home/rstudio/HMMoce_run/compareTracks.r')


myDir <- paste(dataDir, ptt, '/', sep='')
setwd(myDir)
#load(paste(dataDir, ptt, '/', ptt,'_likelihoods.RData', sep=''))
#load(paste(dataDir, ptt, '/', ptt,'_likelihoods.RData', sep=''))
dataDir <- '/home/rstudio/HMMoce_run/data/'
myDir <- paste(dataDir, ptt, '/', sep='')
load(paste(dataDir, 'light_sst_knockout.RData', sep=''))
str(knockout)
#source('~/HMMoce_run/formatTracks.r')

 L.1 <- calc.gpe2(locs, locDates, iniloc = iniloc, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
 L.2 <- calc.sst.par(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
 L.3 <- calc.ohc.par(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
 L.4 <- calc.woa.par(pdt, ptt, woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = F)
 L.5 <- calc.hycom.par(pdt, ptt, hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
 #L.rasters$L.1 <- L.1; 
 #L.rasters$L.2 <- L.2#; L.rasters$L.3 <- L.3; #L.rasters$L.4 <- L.4; 
 #L.rasters$L.5 <- L.5
 resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
 #resamp.idx <- 5
 L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
 
# get idx of layers to keep
#light.idx <- which(dateVec %in% knockout$mako.light)
#sst.idx <- which(dateVec %in% knockout$mako.sst)

# set non-idx layers in L.res to NA
light.knock <- which(dateVec %in% dateVec[-light.idx])
sst.knock <- which(dateVec %in% dateVec[-sst.idx])

for(rr in light.knock){
  L.res[[1]]$L.1[[rr]] <- 0
}

for(rr in sst.knock){
  L.res[[1]]$L.2[[rr]] <- 0
}

if (length(likVec) > 2){
  L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  L.idx <- utils::combn(likVec, 2, simplify=F)
}

#tt <- 13
L.idx <- L.idx[c(1:4, 11:13)]
#L.idx <- L.idx[c(1, 4, 5)]

for (tt in 1:length(L.idx)){
  
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
  
  #bnd <- NA
  i <- 4
  for (bnd in bndVec){
    #for (i in parVec){
      
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
      par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, calcP=T)
      #K1 <- par0$K1; K2 <- par0$K2; 
      P.final <- par0$P.final
      #P.final[1,1] <- .885; P.final[2,2] <- .942
      #P.final[1,2] <- 1-P.final[1,1]; P.final[2,1] <- 1-P.final[2,2]
      # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
      par0 <- makePar(migr.spd=i, grid=g, L.arr=L, calcP=F)
      K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
      
      # RUN THE FILTER STEP
      if(!is.na(bnd)){
        f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
      } else{
        f <- hmm.filter(g, L, K1, K2, P.final, maskL=F)
      }
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K1, K2, L, P.final)
      
      # GET THE MOST PROBABLE TRACK
      tr <- calc.track(s, g, dateVec)
      #plotHMM(s, tr, dateVec, ptt, save.plot = F, behav.pts=T)
      
      # COMPARE HMM, GPE3, SPOT
      setwd(myDir)
      write.table(tr, file=paste(ptt, '_HMM_track.csv', sep=''), sep = ',', col.names = TRUE)
      df <- formatTracks(trackDir = myDir, ptt = ptt, gpeNo = gpeNo, knock=F)
      res <- compareTracks(df)
      res[[4]] <- apply(res$gcd, 2, FUN=function(x) mean(x, na.rm=T))
      res[[5]] <- apply(res$gcd, 2, FUN=function(x) sd(x, na.rm=T))
      res[[6]] <- df
      names(res) <- list('rmse.lon','rmse.lat','gcd','mean.gcd','sd.gcd','tracks')
      #save(res, file=paste(ptt, '_res_knock.RData', sep=''))
      
      # WRITE OUT RESULTS
      outVec <- matrix(c(ptt=paste(ptt,'k_high_f',sep=''), minBounds = bnd, migr.spd = i, rmseLon=res$rmse.lon, rmseLat=res$rmse.lat,
                         gcdMean=res[[4]], gcdSD=res[[5]], paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2]), ncol=30)
      write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
      #save.image(paste(ptt, '_knock_track.RData', sep=''))
      colnames(outVec) <- list('ptt', 'minBnd','migr.spd','rmselon.ti','rmselon.tib','rmselon.kf','rmselon.kfb','rmselon.gpe','rmselon.hmm',
                               'rmselat.ti','rmselat.tib','rmselat.kf','rmselat.kfb','rmselat.gpe','rmselat.hmm',
                               'gcdm.ti','gcdm.tib','gcdm.kf','gcdm.kfb','gcdm.gpe','gcdm.hmm',
                               'gcdsd.ti','gcdsd.tib','gcdsd.kf','gcdsd.kfb','gcdsd.gpe','gcdsd.hmm', 'L.idx', 'migr','resid')
      print(outVec)
      
      
      
    #} # parvec
  } # bnd
} # L.idx tt loop



pdf('check_L_254.pdf', height=18, width=12)
par(mfrow=c(3,2))
for(i in 1:length(dateVec)){
  plot(L.res[[1]]$L.1[[i]]); world(add=T); title(paste(dateVec[i], '-light')); points(spot$lon[i], spot$lat[i])
  plot(L.res[[1]]$L.2[[i]]); world(add=T); title(paste(dateVec[i], '-sst')); points(spot$lon[i], spot$lat[i])
  plot(L.res[[1]]$L.3[[i]]); world(add=T); title(paste(dateVec[i], '-ohc')); points(spot$lon[i], spot$lat[i])
  plot(L.res[[1]]$L.5[[i]]); world(add=T); title(paste(dateVec[i], '-hycom')); points(spot$lon[i], spot$lat[i])
  image.plot(lon,lat,L[i,,]); world(add=T); title(paste(dateVec[i], '-L')); points(spot$lon[i], spot$lat[i])
  image.plot(lon,lat,f$phi[1,i,,]); world(add=T); title(paste(dateVec[i], '-f$phi')); points(spot$lon[i], spot$lat[i])
}

dev.off()


# look at 42

# light 78
  