library(HMMoce)
source('~/Documents/WHOI/RCode/HMMoce_run/formatTracks.r')
source('~/Documents/WHOI/RCode/HMMoce_run/compareTracks.r')

ptts <- c(141254, 141256, 141257, 141259, 121325)

migr <- data.frame(matrix(NA, nrow = 5, ncol = 5))
names(migr) <- c('ptt','m1','m2','m3','m4')
migr$ptt <- ptts
migr$m1 <- c(.5, .5, 1.25, .5, 1.25)
migr$m2 <- c(.75, .75, 1.5, .75, 1.5)
migr$m3 <- c(1, 1, 1.75, 1, 1.75)
migr$m4 <- c(1.25, 1.25, 2, 1.25, 2)

resid <- data.frame(matrix(NA, nrow = 5, ncol = 5))
names(resid) <- c('ptt','r1','r2','r3','r4')
resid$ptt <- ptts
resid$r1 <- c(NA, NA, .5, NA, .5)
resid$r2 <- c(NA, NA, .75, NA, .75)
resid$r3 <- c(NA, NA, 1, NA, 1)
resid$r4 <- c(NA, NA, 1.25, NA, 1.25)

metric.mat <- data.frame(matrix(NA, ncol=7, nrow=80))
names(metric.mat) <- c('ptt','L','migr','resid','mgcd','rmse.lon','rmse.lat')

gpeNo <- c(5, 1, 1, 4, 1)

for (p.idx in 4:length(ptts)){
  ptt <- ptts[p.idx]
  setwd(paste('~/Documents/WHOI/RCode/HMMoce_run/data/', ptt,'/', sep='')) 
  #setwd('C:/RData/HMMoce_run/data/121325/')
  load(paste(ptt,'_hmm_geo.RData',sep=''))
  strt <- p.idx*16-15
  # iterate through L.res1 thru 4
  for (i in 3:4){
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
    
    for (tt in 1:4){
      migr.sel <- migr[which(migr$ptt == ptt), tt+1]
      resid.sel <- resid[which(resid$ptt == ptt), tt+1]
      if(is.na(resid.sel)) resid.sel <- NULL
      
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
      
      par0 <- calc.param(migr.spd = migr.sel, resid.spd = resid.sel, g = g)
      D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
      D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode
      K1 <- gausskern(D1[1], D1[2], muadv = 0)
      K2 <- gausskern(D2[1], D2[2], muadv = 0)
      
      # RUN THE FILTER STEP
      f <- hmm.filter(g, L, K1, K2, P.final)
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K1, K2, P.final)
      
      tr <- calc.track(s, g, dateVec)
      write.table(tr, file=paste(ptt, '_HMM_track.csv', sep=''), sep = ',', col.names = T)
      spot <- read.table(paste(ptt, '_crawl_track.csv', sep=''), sep=',', header=T)
      
      df <- formatTracks(trackDir = getwd(), ptt, gpeNo[p.idx])
      res <- compareTracks(df)
      metrics <- c(mean(res[[3]]$hmm, na.rm=T), res[[1]][6], res[[2]][6]) # MGCD, rmse lon, rmse lat

      base::save.image(paste(ptt,'_hmm_geo.RData', sep=''))
      
      # fill in metrics
      if(is.null(resid.sel)) resid.sel <- NA
      metric.mat[strt.tt,] <- c(ptt, i, migr.sel, resid.sel, metrics)
      print(metric.mat[strt.tt,])
      #then add to strt.tt
      strt.tt <- strt.tt + 1
    }
    
  }
  
}
