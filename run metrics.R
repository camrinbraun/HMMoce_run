library(HMMoce)
metric.mat <- data.frame(matrix(NA, ncol=10, nrow=40))

gpeNo <- c(7, 5, 4, 7)
ptts <- c(141254, 141256, 141257, 141259)

for (p.idx in 2:length(ptts)){
  ptt <- ptts[p.idx]
  setwd(paste('~/Documents/WHOI/RCode/HMMoce_run/data/', ptt,'/', sep='')) 
  #setwd('C:/RData/HMMoce_run/data/121325/')

  # load knockout
  load(paste(ptt,'_hmm_geo_knockout.RData', sep=''))
  
  # change L names
  L.sst.k <- L.sst.par; L.light.k <- L.light
  
  # load full RData
  load(paste(ptt,'_hmm_geo.RData',sep=''))
  strt <- p.idx*16-15
  
  # create a list of all the desired input likelihood rasters
  L.rasters1 <- list(L.sst = L.sst.k, L.light = L.light.k)
  L.rasters2 <- list(L.sst = L.sst.k, L.light = L.light.k, L.ohc = L.ohc)
  L.rasters3 <- list(L.sst = L.sst.k, L.light = L.light.k, L.prof = L.hycom.par)
  L.rasters4 <- list(L.sst = L.sst.k, L.light = L.light.k, L.prof = L.qwoa.par)
  
  # L.sst is the resolution/extent we're sampling everything TO
  L.res1 <- resample.grid.par(L.rasters1, L.rasters1$L.sst)
  L.res2 <- resample.grid.par(L.rasters2, L.rasters2$L.sst)
  L.res3 <- resample.grid.par(L.rasters3, L.rasters3$L.sst)
  L.res4 <- resample.grid.par(L.rasters4, L.rasters4$L.sst)
  
  gpeNo <- c(7, 5, 4, 7)
  metric.mat <- data.frame(matrix(NA, ncol=10, nrow=40))
  ptts <- c(141254, 141256, 141257, 141259)
  
  # iterate through L.res1 thru 4
  for (i in 1:4){
    if(i==1){
      L.mle.res <- L.res1$L.mle.res
      g <- L.res1$g; lon <- g$lon[1,]; lat <- g$lat[,1]
      g.mle <- L.res1$g.mle
      
      L <- make.L(L1 = L.res1[[1]]$L.sst,
                  L2 = L.res1[[1]]$L.light,
                  L.mle.res = L.mle.res, dateVec = dateVec,
                  locs.grid = locs.grid, iniloc = iniloc)
      
      L.mle <- L$L.mle; L <- L$L
    } else if(i == 2){
      L.mle.res <- L.res2$L.mle.res
      g <- L.res2$g; lon <- g$lon[1,]; lat <- g$lat[,1]
      g.mle <- L.res2$g.mle
      
      L <- make.L(L1 = L.res2[[1]]$L.sst,
                  L2 = L.res2[[1]]$L.light,
                  L3 = L.res2[[1]]$L.ohc,
                  L.mle.res = L.mle.res, dateVec = dateVec,
                  locs.grid = locs.grid, iniloc = iniloc)
      
      L.mle <- L$L.mle; L <- L$L
    } else if(i == 3){
      L.mle.res <- L.res3$L.mle.res
      g <- L.res3$g; lon <- g$lon[1,]; lat <- g$lat[,1]
      g.mle <- L.res3$g.mle
      
      L <- make.L(L1 = L.res3[[1]]$L.sst,
                  L2 = L.res3[[1]]$L.light,
                  L3 = L.res3[[1]]$L.prof,
                  L.mle.res = L.mle.res, dateVec = dateVec,
                  locs.grid = locs.grid, iniloc = iniloc)
      
      L.mle <- L$L.mle; L <- L$L
    } else if(i == 4){
      L.mle.res <- L.res4$L.mle.res
      g <- L.res4$g; lon <- g$lon[1,]; lat <- g$lat[,1]
      g.mle <- L.res4$g.mle
      
      L <- make.L(L1 = L.res4[[1]]$L.sst,
                  L2 = L.res4[[1]]$L.light,
                  L3 = L.res4[[1]]$L.prof,
                  L.mle.res = L.mle.res, dateVec = dateVec,
                  locs.grid = locs.grid, iniloc = iniloc)
      
      L.mle <- L$L.mle; L <- L$L
    }
    
    if (max(L[1,,], na.rm=T) > 1){
      print('L[1,,] > 1')
    } else{
      stop('Error. L tag and pop not greater than 1.')
      print('Error. L tag and pop not greater than 1.')
      break
      
    }
    
    if(exists('strt.tt')){
      strt.tt <- strt.tt+1
    } else{
      strt.tt <- strt+i-1
    }
    
    #for (tt in 1:4){
      migr.sel <- 2
      resid.sel <- NULL
      
      par0 <- calc.param(migr.spd = migr.sel, resid.spd = resid.sel, g = g.mle)
      if(par0[[1]] <= .5) par0[[1]] <- .51
      if(par0[[3]] <= .5) par0[[3]] <- .51
      D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
      D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode
      
      # GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
      K1 <- gausskern(D1[1], D1[2], muadv = 0)
      K2 <- gausskern(D2[1], D2[2], muadv = 0)
      
      # MAKE A GUESS AT STATE SWITCHING PROBABILITY
      p <- c(0.7, 0.8)
      
      # RUN EXPECTATION-MAXIMIZATION ROUTINE FOR MATRIX, P (STATE SWITCH PROBABILITY)
      P.init <- matrix(c(p[1], 1 - p[1], 1 - p[2], p[2]), 2, 2, byrow = TRUE)
      P.final <- expmax(P.init, g = g.mle, L = L.mle, K1, K2, save = T)
      save.p <- P.final[[2]]; P.final <- P.final[[1]]
      
      par0 <- calc.param(migr.spd = 4, resid.spd = resid.sel, g = g)
      D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
      D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode
      K1 <- gausskern(D1[1], D1[2], muadv = 0)
      K2 <- gausskern(D2[1], D2[2], muadv = 0)
      
      # RUN THE FILTER STEP
      f <- hmm.filter(g, L, K1, K2, P.final)
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K1, K2, P.final)
      
      tr <- calc.track(s, g, dateVec)
      write.table(tr, file=paste(ptt, '_HMM_track.csv', sep=''), sep = ',')#, col.names = T)
      spot <- read.table(paste(ptt, '_crawl_track.csv', sep=''), sep=',', header=T)
      
      df <- formatTracks(trackDir = getwd(), ptt, gpeNo[p.idx])
      res <- compareTracks(df)
      metrics <- c(mean(res[[3]]$hmm, na.rm=T), res[[1]][6], res[[2]][6],
                   mean(res[[3]]$gpe, na.rm=T), res[[1]][5], res[[2]][5]) # MGCD, rmse lon, rmse lat
      
      #base::save.image(paste(ptt,'_hmm_geo_k1.RData', sep=''))
      
      # fill in metrics
      if(is.null(resid.sel)) resid.sel <- NA
      metric.mat[strt.tt,] <- c(ptt, i, migr.sel, resid.sel, metrics)
      print(metric.mat[strt.tt,])
      #then add to strt.tt
      strt.tt <- strt.tt + 1
    }
    
  }
  
}

plot(df$spot.lon, df$spot.lat, type='l')
lines(df$gpe.lon, df$gpe.lat, col='blue')
lines(df$hmm.lon, df$hmm.lat, col='green')

pdf('try_plots.pdf', height=12, width=8)
par(mfrow=c(2,1))
for(i in 1:181){
  image.plot(lon,lat, L[i,,]); points(spot$lon[i], spot$lat[i], col='white')
  image.plot(lon,lat,f$phi[1,i,,]); points(spot$lon[i], spot$lat[i], col='white')
  
}
dev.off()