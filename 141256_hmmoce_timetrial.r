# RUN BLUE 259 VIA HMMoce
library(HMMoce)

# SETWD
setwd('~/ebs/Data/BlueSharks/141256/') 
#setwd('C:/RData/HMMoce_run/data/141256/')
load('141256_likelihoods.RData')


#----------------------------------------------------------------------------------#
# CALC LIKELIHOODS
#----------------------------------------------------------------------------------#

# GENERATE LIGHT LIKELIHOOD
# SRSS METHOD
#L.light <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec)
# OR
# GPE2 METHOD
t001 <- Sys.time()

L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = F, gpeOnly = TRUE)

# GENERATE DAILY SST LIKELIHOODS
#t0 <- Sys.time()
sst.dir <- '/home/rstudio/ebs/EnvData/sst/BlueSharks/141256/'
L.2 <- calc.sst.par(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec, sens.err = 2.5)
#L.2 <- calc.sst.par(tag.sst, filename=fname, sst.dir = sst.dir.small, dateVec = dateVec, sens.err = 1)

#t1 <- Sys.time()
#L.sst.par <- calc.sst.par(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec, sens.err = 2.5)
#t2 <- Sys.time()

#-------
# GENERATE DAILY OCEAN HEAT CONTENT (OHC) LIKELIHOODS
#pdt.try <- pdt[-which(pdt$Date == '2016-01-03 00:00:00'),]
t0 <- Sys.time()
hycom.dir <- '/home/rstudio/ebs/EnvData/hycom/'
L.3 <- calc.ohc.par(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = FALSE)
#L.3 <- calc.ohc.par(pdt, filename=fname, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F, ncores=parallel::detectCores()-2)
t1 <- Sys.time()
#L.ohc.par <- calc.ohc.par(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = FALSE)
t2 <- Sys.time()

# WOA DATA
#woa.dir <- getwd()#paste('my_woa_dir')
#get.env(type = 'woa', resol = 'one')

#t3 <- Sys.time()
#L.woa.par <- calc.woa.par(pdt, ptt, woa.data = woa.one, focalDim = 3, dateVec = dateVec, use.se = F)
#t4 <- Sys.time()
#setwd(myDir)
#base::save.image('141256_hmm_geo.RData')

#t1<-Sys.time()
#L.qwoa.try <- calc.profile(pdt, ptt, dat = woa.quarter$watertemp, lat = woa.quarter$lat, lon = woa.quarter$lon, dateVec = dateVec, envType = 'woa')
#t5 <- Sys.time()
#L.qwoa.par <- calc.woa.par(pdt.try, ptt, woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = F)
#t6 <- Sys.time()

t7<-Sys.time()
L.5 <- calc.hycom.par(pdt, ptt, hycom.dir, focalDim = 9, dateVec = dateVec, use.se = F)
t8<-Sys.time()
#L.hycom <- calc.profile(pdt, ptt, hycom.dir = hycom.dir, dateVec = dateVec, envType = 'hycom')
#t3<-Sys.time()
setwd(myDir)
#base::save.image('141256_hmm_geo_v5.RData')

t002 <- Sys.time()

gc(); closeAllConnections()

#----------------------------------------------------------------------------------#
# LIST, RESAMPLE, SAVE
#----------------------------------------------------------------------------------#
#for (iii in likVec){
#  if(!exists(paste('L.',iii,sep=''))){
##    statusVec <- c(statusVec, paste('Error: L.', iii, ' does not exist. Cannot combine required likelihoods.',sep=''))
#    stop(paste('Error: L.', iii, ' does not exist. Cannot combine required likelihoods.',sep=''))
#  } 
#}

# eventually can be used to fill any empty L in L.res
#rvec <- list()
#for (iii in 1:5){
#  rvec[[iii]] <- exists(paste('L.',iii,sep=''))
#}
t003 <- Sys.time()
rm(L.4); rm(L.idx); rm(L.rasters); rm(L.res)
L.rasters <- mget(ls(pattern = 'L\\.'))
#resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
L.res <- resample.grid(L.rasters, raster::aggregate(L.rasters[[1]], 4), mle.res=1)
#L.res <- resample.grid(L.rasters, aggregate(L.1, 3), mle.res=1)
#L.res <- resample.grid(L.rasters, aggregate(L.1, 2), mle.res=1)
#L.res <- resample.grid(L.rasters, L.1, mle.res=1)
#L.res <- resample.grid(L.rasters, L.3, mle.res=1)


# Figure out appropriate L combinations
if (length(likVec) > 2){
  L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  L.idx <- utils::combn(likVec, 2, simplify=F)
}

rm(list=c('L.1','L.2','L.3','L.4','L.5', 'woa.quarter'))

# save workspace image to s3 as checkpoint
#setwd(myDir); base::save.image('check2.rda')
#aws.s3::s3save_image(bucket=paste(bucketDir, '/', ptt, sep=''), object='check2.rda')
#if(file.exists('check2.rda')) file.remove('check1.rda')

#================
## END STEP 2
#================
#enterAt <- 3
#} 
##rm(list=ls())
#}

#if (enterAt == 3){
#dataDir <- 'C:/RData/HMMoce_run/data/' 
#myDir <- paste(dataDir, ptt, '/', sep='')
#  setwd(myDir)#; load('check2.rda')
#t003 <- Sys.time()
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
    #setwd('~/HMMoce'); devtools::load_all()
    #setwd(myDir)
    library(HMMoce)
    #for (tt in run.idx){
  #bnd <- NA
 # i=4
#  
  for (bnd in bndVec){
      for (i in parVec){
        
        runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i,sep='')
        
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
        save(res, file=paste(runName, '-HMMoce_res.rda', sep=''))
        #save.image(file=paste(ptt, '-HMMoce.RData', sep=''))
        #aws.s3::s3save(res, bucket=paste(bucketDir, '/', ptt, sep=''), object=paste(runName, '-HMMoce_res.rda', sep=''))
        #source('~/HMMoce/R/hmm.diagnose.r')
        #hmm.diagnose(res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt, plot=T)
        
        #outVec <- outVec
        
      } # parVec loop
    } # bndVec loop
  } # L.idx loop
    
  
  parallel::stopCluster(cl)
  closeAllConnections()
  


t004 <- Sys.time()
      
#-------
# GENERATE DAILY PROFILE LIKELIHOODS
#L.prof.woa <- calc.profile(pdt, dat = woa, lat = lat, lon = lon, dateVec = dateVec, envType = 'woa')

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#
# create a list of all the desired input likelihood rasters
L.rasters1 <- list(L.sst = L.sst, L.light = L.light)
L.rasters2 <- list(L.sst = L.sst, L.light = L.light, L.ohc = L.ohc)
L.rasters3 <- list(L.sst = L.sst, L.ohc = L.ohc)
#L.rasters3 <- list(L.sst = L.sst, L.light = L.light, L.prof = L.hycom.par)
#L.rasters4 <- list(L.sst = L.sst, L.light = L.light, L.prof = L.qwoa.par)

# L.sst is the resolution/extent we're sampling everything TO
L.res1 <- resample.grid.par(L.rasters1, L.rasters1$L.sst)
L.res2 <- resample.grid.par(L.rasters2, L.rasters2$L.sst)
L.res3 <- resample.grid.par(L.rasters3, L.rasters3$L.sst)
#L.res4 <- resample.grid.par(L.rasters4, L.rasters4$L.sst)

#setwd(myDir)
#base::save.image('141256_hmm_geo.RData')

# pull some other helpful variables from the resample.grid() output for later use
L.mle.res <- L.res1$L.mle.res
g <- L.res1$g; lon <- g$lon[1,]; lat <- g$lat[,1]
g.mle <- L.res$g.mle

#----------------------------------------------------------------------------------#
# LOAD AND FORMAT DATAFRAME OF KNOWN LOCATIONS, IF ANY
#----------------------------------------------------------------------------------#

#colnames(known.locs) <- list('date','lat','lon')
#   where 'date' is from as.Date(known.locs$date)

#----------------------------------------------------------------------------------#
# COMBINE LIKELIHOOD MATRICES
#----------------------------------------------------------------------------------#
# this example just uses L.sst and L.light. You can list up to three (L1, L2, L3 inputs).
L <- make.L.mod(L1 = L.res1[[1]]$L.sst,
            L2 = L.res1[[1]]$L.light,
            #L3 = L.res2[[1]]$L.ohc,
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc, bathy=bathy,
            pdt = pdt)

L.mle <- L$L.mle; L <- L$L

#----------------------------------------------------------------------------------#
# FIGURE OUT MOVEMENT PARAMETERS
#----------------------------------------------------------------------------------#

# PROVIDE FIXED KERNEL PARAMETERS
par0 <- calc.param2(migr.spd = 2, g = g)
#par0$resid <- par0$migr; par0$sig2 <- par0$sig1/10
#par0 <- c(8.908,10.27,1.152,0.0472)
#par0 <- c(6.2, 8, 6.2, 0.05)
D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode

# GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
K1 <- gausskern(D1[1], D1[2], muadv = 0)
K2 <- gausskern(D2[1], D2[2], muadv = 0)

# MAKE A GUESS AT STATE SWITCHING PROBABILITY
p <- c(0.7, 0.8)

# RUN EXPECTATION-MAXIMIZATION ROUTINE FOR MATRIX, P (STATE SWITCH PROBABILITY)
P.init <- matrix(c(p[1], 1 - p[1], 1 - p[2], p[2]), 2, 2, byrow = TRUE)
#P.final <- expmax(P.init, g = g.mle, L = L.mle, K1, K2, save = T)
#save.p <- P.final[[2]]; P.final <- P.final[[1]]
P.final <- P.init
#----------------------------------------------------------------------------------#
#par0 <- calc.param(migr.spd = 2, g = g)
#D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
#D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode
#K1 <- gausskern(D1[1], D1[2], muadv = 0)
#K2 <- gausskern(D2[1], D2[2], muadv = 0)

# RUN THE FILTER STEP
f <- hmm.filter.ext(g, L, K1, K2, maskL=T, par0=par0, P.final)
f.nom <- hmm.filter.ext(g, L, K1, K2, maskL=F, P.final)
f.old <- hmm.filter(g, L, K1, K2, P.final)

#pdf('check_filter_oldfilter.pdf', height=8, width=12)
#par(mfrow=c(1,2))
#lon <- g$lon[1,]; lat <- g$lat[,1]
#for(i in 1:length(dateVec)){
#  image.plot(lon,lat,L[i,,]); points(spot$lon[i], spot$lat[i]);world(add=T); title(paste('L-',dateVec[i],sep=''))
#  image.plot(lon,lat,f$pred[1,i,,]);points(spot$lon[i], spot$lat[i]); world(add=T); title('f$pred[1,,]')
#}
#dev.off()
# plot if you want to see confidence limits
#res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
#fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P.final)
s.old = hmm.smoother(f.old, K1, K2, P.final)

# plot if you want to see confidence limits
#sres = apply(s[1,,,], 2:3, sum, na.rm=T)
#fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#

tr <- calc.track(s, g, dateVec)

plotHMM(s, tr, dateVec, ptt, save.plot = F)

plot(tr$lon,tr$lat,type='l');world(add=T);lines(spot$lon,spot$lat,col='red')

#setwd(myDir)
write.table(tr, file=paste(ptt, '_HMM_track.csv', sep=''), sep = ',', col.names = T)
base::save.image(paste(ptt,'_hmm_geo_v5.RData', sep=''))

#=======================================================================================#
## END
#=======================================================================================#

spot <- read.table('141256_crawl_track.csv', sep=',', header = T)

pdf('check_141256_L.pdf',width=12,height=16)
par(mfrow=c(2,2))
for (i in 1:length(dateVec)){
  plot(L.res2[[1]]$L.light[[i]]); world(add=T); points(spot$lon[i], spot$lat[i], pch=21); title(paste('Light-', dateVec[i],sep=''))
  plot(L.res2[[1]]$L.sst[[i]]); world(add=T); points(spot$lon[i], spot$lat[i], pch=21); title(paste('SST-', dateVec[i],sep=''))
  plot(L.res2[[1]]$L.ohc[[i]]); world(add=T); points(spot$lon[i], spot$lat[i], pch=21); title(paste('OHC-', dateVec[i],sep=''))
  image.plot(lon, lat, L[i,,]); world(add=T, fill=T); points(spot$lon[i], spot$lat[i], col='white', pch=21); title(paste('L-', dateVec[i],sep=''))
  #print(i)
  #image.plot(lon, lat, f$pred[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  #image.plot(lon, lat, s[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  #image.plot(lon, lat, s[2,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
}
dev.off()

pdf('check_filter_141256_old.pdf',width=16,height=12)
par(mfrow=c(2,2))
for (i in 1:length(dateVec)){
  #plot(L.res2[[1]]$L.light[[i]]); world(add=T); #points(spot$lon[i], spot$lat[i], pch=21); title(paste('Light-', dateVec[i],sep=''))
  #plot(L.res2[[1]]$L.sst[[i]]); world(add=T); #points(spot$lon[i], spot$lat[i], pch=21); title(paste('SST-', dateVec[i],sep=''))
  #plot(L.res2[[1]]$L.ohc[[i]]); world(add=T); #points(spot$lon[i], spot$lat[i], pch=21); title(paste('OHC-', dateVec[i],sep=''))
  image.plot(lon, lat, L[i,,]); world(add=T, fill=T); points(spot$lon[i], spot$lat[i], col='white', pch=21); title(paste('L-', dateVec[i],sep=''))
  #print(i)
  image.plot(lon, lat, f$pred[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i]); title('filter$pred')
  image.plot(lon, lat, s[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  #image.plot(lon, lat, s[2,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  image.plot(lon, lat, f$phi[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i]); title('filter$phi')
  
}
dev.off()

L.ohc.mod <- L.ohc / max(L.ohc, na.rm=T) * 2
plot(L.ohc.mod[[123]])
plot(L.ohc.mod[[121]])
L.ohc.mod[L.ohc.mod == 0] <- 1
L.ohc.mod[is.na(L.ohc.mod)] <- 0
try.r <- prod(L.ohc.mod[[123]], L.ohc.mod[[121]], na.rm=F)
try.r[try.r == 0] <- NA
try.r[try.r == 1] <- 0

plot(try.r)


sumIdx <- which(raster::cellStats(L.light, sum, na.rm = T) != 0)
for (i in sumIdx){
  L.light[[i]] <- L.light[[i]] / raster::cellStats(L.light[[i]], max, na.rm = T)
}


## need an example to illustrate the problem(s) with filter:

# start with a known location and simple migratory-like diffusion:
image.plot(lon,lat,L[1,,])
world(add=T)

K1 <- gausskern(6, 6, muadv = 0)
image.plot(K1)
K1 <- imager::as.cimg(K1)


# then apply some of the filter logic:
# convolve the known location to get t+1
#p1 = imager::as.cimg(t(phi[1, t-1,,]))
p1 = imager::as.cimg(t(L[1,,]))
q1 = imager::convolve(p1, K1)
q1 = t(as.matrix(q1))

# multiply by data likelihood, L[2,,] has good data so let's cause the problem with different day
t=2
post1 <- q1 * L[t,,]

# so the 

