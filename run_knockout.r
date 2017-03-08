pttList <- c(141259, 141257, 141256, 141254)
load('light_sst_knockout.RData')

for (tt in 1:length(pttList)){
  
  ptt <- pttList[tt]
  
  dataDir <- paste('~/Documents/WHOI/RCode/HMMoce_run/data/', ptt, sep='')
  envDir <- paste('~/Documents/WHOI/RCode/HMMoce_run/env_data/', ptt, sep='')
  
  # TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
  inilocList <- list(data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 
                                         10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T)),
                     data.frame(matrix(c(15, 10, 2015, 41.637, -69.706, 
                                         12, 4, 2016, 37.751, -71.649), nrow = 2, ncol = 5, byrow = T)),
                     data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                                         24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T)),
                     data.frame(matrix(c(21, 10, 2015, 41.597,	-69.445, 
                                         5, 2, 2016, 38.589, -54.874), nrow = 2, ncol = 5, byrow = T)))
  for (i in 1:4){
    names(inilocList[[i]]) <- list('day','month','year','lat','lon')
  }
  
  iniloc <- inilocList[[tt]]
  
  sp.limList <- list(list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50),
                     list(lonmin = -85, lonmax = -50, latmin = 20, latmax = 50),
                     list(lonmin = -95, lonmax = -52, latmin = 10, latmax = 55),
                     list(lonmin = -85, lonmax = -45, latmin = 30, latmax = 50))
  
  sp.lim <- sp.limList[[tt]]
  
  runHMM(ptt, iniloc, dataDir, envDir, sp.lim)
  
}

runHMM <- function(ptt, iniloc, dataDir, envDir, sp.lim){
  
  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day')) 
  
  # READ IN DATA FROM WC FILES
  #myDir <- '~/Documents/WHOI/RCode/HMMoce/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
  myDir <- dataDir
  # sst data
  tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
  sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data
  tag.sst$dts <- as.Date(as.POSIXct(tag.sst$Date, format=HMMoce:::findDateFormat(tag.sst$Date)))
  tag.sst <- tag.sst[which(tag.sst$dts %in% knockout$blue259.sst),]
  sst.udates <- unique(tag.sst$dts)
  tag.sst <- tag.sst[,c(1:11)]
  
  # depth-temp profile data
  pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
  pdt.udates <- pdt$udates; pdt <- pdt$data
  
  # light data
  light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
  light.udates <- light$udates; light <- light$data
  light$dts <- as.Date(as.POSIXct(light$Day, format='%d-%b-%y'))
  light <- light[which(light$dts %in% knockout$blue259.light),]
  light.udates <- unique(light$dts)
  light <- light[,c(1:51)]
  
  # OPTIONAL: light data as output from GPE2, different filtering algorithm seems to work better for light likelihood generation
  #locs <- read.table(paste(myDir, '/', ptt, '-Locations-GPE2.csv', sep = ''), sep = ',', header = T, blank.lines.skip = F)
  #locDates <- as.Date(locs$Date)
  
  #----------------------------------------------------------------------------------#
  # FURTHER PREPARATION
  # Set spatial limits and download env data
  #----------------------------------------------------------------------------------#
  
  # SET SPATIAL LIMITS, IF DESIRED
  #sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)
  
  if (exists('sp.lim')){
    locs.grid <- setup.locs.grid(sp.lim)
  } else{
    locs.grid <- setup.locs.grid(locs)
    sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                   latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
  }
  
  # IF YOU NEED TO DOWNLOAD SST DATA
  sst.dir <- paste(envDir, '/sst/',sep='')
  #sst.dir <- '~/Documents/WHOI/RCode/HMMoce_run/env_data/141259/sst/'
  #sst.dir <- 'C:/RData/HMMoce_run/env_data/141259/sst/'
  #get.env(sst.udates, type = 'sst', spatLim = sp.lim, save.dir = sst.dir)
  
  # HYCOM DATA
  hycom.dir <- paste(envDir, '/hycom/', sep='')
  #hycom.dir <- paste('C:/RData/HMMoce_run/env_data/141259/hycom/')
  #get.env(pdt.udates, type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
  
  # AND/OR WOA DATA
  woa.dir <- '~/Documents/WHOI/RCode/HMMoce_run/env_data/woa/'
  #woa.dir <- 'C:/RData/HMMoce_run/env_data/woa/'
  #get.env(type = 'woa', resol = 'quarter')
  # then load the downloaded rda file
  load(paste(woa.dir,'woa.quarter.rda',sep=''))
  #str(woa.quarter)
  #List of 4
  #$ watertemp: num [1:44, 1:46, 1:57, 1:12] 26.5 26.5 26.4 26.3 26.2 ...
  #$ lon      : num [1:44(1d)] -95.5 -94.5 -93.5 -92.5 -91.5 -90.5 -89.5 -88.5 -87.5 -86.5 ...
  #$ lat      : num [1:46(1d)] 9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5 ...
  #$ depth    : num [1:57(1d)] 0 5 10 15 20 25 30 35 40 45 ...
  
  #----------------------------------------------------------------------------------#
  # CALC LIKELIHOODS
  #----------------------------------------------------------------------------------#
  
  # GENERATE LIGHT LIKELIHOOD
  # SRSS METHOD
  L.light <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec)
  # OR
  # GPE2 METHOD
  #L.locs <- calc.gpe2(locs, locDates, iniloc = iniloc, locs.grid = locs.grid, dateVec = dateVec, errEll = TRUE, gpeOnly = TRUE)
  
  
  # GENERATE DAILY SST LIKELIHOODS
  #L.sst <- calc.sst(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec)
  #t1 <- Sys.time()
  L.sst.par <- calc.sst.par(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec)
  #t2 <- Sys.time()
  
  #-------
  # GENERATE DAILY OCEAN HEAT CONTENT (OHC) LIKELIHOODS
  #t1 <- Sys.time()
  L.ohc <- calc.ohc.par(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '')
  #t2 <- Sys.time()
  
  # WOA DATA
  #t3 <- Sys.time()
  #L.woa.par <- calc.woa.par(pdt, ptt, woa.data = woa.one, focalDim = 3, dateVec = dateVec)
  #t4 <- Sys.time()
  
  #t1<-Sys.time()
  #L.qwoa.try <- calc.profile(pdt, ptt, dat = woa.quarter$watertemp, lat = woa.quarter$lat, lon = woa.quarter$lon, dateVec = dateVec, envType = 'woa')
  #t5<-Sys.time()
  L.qwoa.par <- calc.woa.par(pdt, ptt, woa.data = woa.quarter, focalDim = 9, dateVec = dateVec)
  #t6<-Sys.time()
  
  #t7<-Sys.time()
  L.hycom.par <- calc.hycom.par(pdt, ptt, hycom.dir, focalDim = 9, dateVec = dateVec)
  #t8<-Sys.time()
  #L.hycom <- calc.profile(pdt, ptt, hycom.dir = hycom.dir, dateVec = dateVec, envType = 'hycom')
  #t3<-Sys.time()
  
  #-------
  # GENERATE DAILY PROFILE LIKELIHOODS
  #L.prof.woa <- calc.profile(pdt, dat = woa, lat = lat, lon = lon, dateVec = dateVec, envType = 'woa')
  
  #----------------------------------------------------------------------------------#
  # SETUP A COMMON GRID
  #----------------------------------------------------------------------------------#
  # create a list of all the desired input likelihood rasters
  L.rasters1 <- list(L.sst = L.sst.par, L.locs = L.locs)
  L.rasters2 <- list(L.sst = L.sst.par, L.locs = L.locs, L.ohc = L.ohc)
  #L.rasters3 <- list(L.sst = L.sst.par, L.locs = L.locs, L.prof = L.hycom.par)
  #L.rasters4 <- list(L.sst = L.sst.par, L.locs = L.locs, L.prof = L.qwoa.par)
  
  # L.sst is the resolution/extent we're sampling everything TO
  L.res1 <- resample.grid.par(L.rasters1, L.rasters1$L.sst)
  L.res2 <- resample.grid.par(L.rasters2, L.rasters2$L.sst)
  #L.res3 <- resample.grid.par(L.rasters3, L.rasters3$L.prof)
  #L.res4 <- resample.grid.par(L.rasters4, L.rasters4$L.sst)
  
  # pull some other helpful variables from the resample.grid() output for later use
  L.mle.res <- L.res2$L.mle.res
  g <- L.res2$g; lon <- g$lon[1,]; lat <- g$lat[,1]
  g.mle <- L.res2$g.mle
  
  #----------------------------------------------------------------------------------#
  # LOAD AND FORMAT DATAFRAME OF KNOWN LOCATIONS, IF ANY
  #----------------------------------------------------------------------------------#
  
  #colnames(known.locs) <- list('date','lat','lon')
  #   where 'date' is from as.Date(known.locs$date)
  
  #----------------------------------------------------------------------------------#
  # COMBINE LIKELIHOOD MATRICES
  #----------------------------------------------------------------------------------#
  # this example just uses L.sst and L.light. You can list up to three (L1, L2, L3 inputs).
  L <- make.L(L1 = L.res2[[1]]$L.sst,
              L2 = L.res2[[1]]$L.locs,
              L3 = L.res2[[1]]$L.ohc,
              L.mle.res = L.mle.res, dateVec = dateVec,
              locs.grid = locs.grid, iniloc = iniloc)
  
  L.mle <- L$L.mle; L <- L$L
  
  #----------------------------------------------------------------------------------#
  # FIGURE OUT MOVEMENT PARAMETERS
  #----------------------------------------------------------------------------------#
  # PROVIDE FIXED KERNEL PARAMETERS
  par0 <- calc.param(migr.spd = 2, g = g.mle)
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
  P.init <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)
  P.final <- expmax(P.init, g = g.mle, L = L.mle, K1, K2)
  #save.p <- P.final[[2]]; P.final <- P.final[[1]]
  
  #----------------------------------------------------------------------------------#
  # re-run parameter setup to adjust for using g, rather than g.mle as in P = expmax() above
  par0 <- calc.param(migr.spd = 2, g = g)
  D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
  D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode
  K1 <- gausskern(D1[1], D1[2], muadv = 0)
  K2 <- gausskern(D2[1], D2[2], muadv = 0)
  
  # RUN THE FILTER STEP
  f <- hmm.filter(g, L, K1, K2, P.final)
  
  # plot if you want to see confidence limits
  #res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
  #fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))
  
  #----------------------------------------------------------------------------------#
  # RUN THE SMOOTHING STEP
  s = hmm.smoother(f, K1, K2, P.final)
  
  # plot if you want to see confidence limits
  #sres = apply(s[1,,,], 2:3, sum, na.rm=T)
  #fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))
  
  #----------------------------------------------------------------------------------#
  # GET THE MOST PROBABLE TRACK
  #----------------------------------------------------------------------------------#
  
  tr <- calc.track(s, g, dateVec)
  
  plotHMM(s, tr, dateVec, ptt, save.plot = F)
  #setwd(myDir)
  write.table(tr, file=paste(ptt, '_HMM_track_knockout.csv', sep=''), sep = ',', col.names = T)
  base::save.image(paste(ptt,'_hmm_geo_knockout.RData', sep=''))
  
  #=======================================================================================#
  ## END
  #=======================================================================================#
  
  
}

