
# RUN BLUE 259 VIA HMMoce
library(HMMoce)

# SETWD
setwd('~/Documents/WHOI/RCode/HMMoce_run/data/141259/') 
#setwd('C:/RData/HMMoce_run/data/141259/')
load('141259_hmm_geo.RData')

# READ IN TAG DATA
ptt <- 141259

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 
                              10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
#myDir <- '~/Documents/WHOI/RCode/HMMoce/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
myDir <- getwd()
# sst data
tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

# depth-temp profile data
pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
pdt.udates <- pdt$udates; pdt <- pdt$data

# light data
light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
light.udates <- light$udates; light <- light$data

# OPTIONAL: light data as output from GPE2, different filtering algorithm seems to work better for light likelihood generation
locs <- read.table(paste(ptt, '-Locations-GPE2.csv', sep = ''), sep = ',', header = T, blank.lines.skip = F)
locDates <- as.Date(locs$Date)

#----------------------------------------------------------------------------------#
# FURTHER PREPARATION
# Set spatial limits and download env data
#----------------------------------------------------------------------------------#

# SET SPATIAL LIMITS, IF DESIRED
sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)

if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}

# IF YOU NEED TO DOWNLOAD SST DATA
sst.dir <- '~/Documents/WHOI/RCode/HMMoce_run/env_data/141259/sst/'
#sst.dir <- 'C:/RData/HMMoce_run/env_data/141259/sst/'
#get.env(sst.udates, ptt = ptt, type = 'sst', spatLim = sp.lim, save.dir = sst.dir)

# HYCOM DATA
hycom.dir <- '~/Documents/WHOI/RCode/HMMoce_run/env_data/141259/hycom/'
#hycom.dir <- paste('C:/RData/HMMoce_run/env_data/141259/hycom/')
#get.env(pdt.udates[2:length(pdt.udates)], ptt = ptt, type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

# AND/OR WOA DATA
woa.dir <- '~/Documents/WHOI/RCode/HMMoce_run/env_data/woa/'
#woa.dir <- 'C:/RData/HMMoce_run/env_data/woa/'
#get.env(type = 'woa', resol = 'quarter')
# then load the downloaded rda file
load(paste(woa.dir,'woa.quarter.rda',sep=''))
str(woa.quarter)
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
#L.light <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec)
# OR
# GPE2 METHOD
L.light <- calc.gpe2(locs, locDates, iniloc = iniloc, locs.grid = locs.grid, dateVec = dateVec, errEll = F, gpeOnly = TRUE)


# GENERATE DAILY SST LIKELIHOODS
#L.sst <- calc.sst(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec, sens.err = 1)
t1 <- Sys.time()
L.sst.par <- calc.sst.par(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec, sens.err = 2.5)
t2 <- Sys.time()

#-------
# GENERATE DAILY OCEAN HEAT CONTENT (OHC) LIKELIHOODS
#t3 <- Sys.time()
#L.ohc <- calc.ohc(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = FALSE)
t3 <- Sys.time()
L.ohc.par <- calc.ohc.par(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = FALSE)
t4 <- Sys.time()

# WOA DATA
#woa.dir <- getwd()#paste('my_woa_dir')
#get.env(type = 'woa', resol = 'one')

#t3 <- Sys.time()
#L.woa.par <- calc.woa.par(pdt, ptt, woa.data = woa.one, focalDim = 3, dateVec = dateVec)
#t4 <- Sys.time()

#t1<-Sys.time()
#L.qwoa.try <- calc.profile(pdt, ptt, dat = woa.quarter$watertemp, lat = woa.quarter$lat, lon = woa.quarter$lon, dateVec = dateVec, envType = 'woa')
t5<-Sys.time()
L.qwoa.par <- calc.woa.par(pdt, ptt, woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = T)
t6<-Sys.time()

t7<-Sys.time()
L.hycom.par <- calc.hycom.par(pdt, ptt, hycom.dir, focalDim = 9, dateVec = dateVec, use.se = F)
t8<-Sys.time()
#L.hycom <- calc.profile(pdt, ptt, hycom.dir = hycom.dir, dateVec = dateVec, envType = 'hycom')
#t3<-Sys.time()

#-------
# GENERATE DAILY PROFILE LIKELIHOODS
#L.prof.woa <- calc.profile(pdt, dat = woa, lat = lat, lon = lon, dateVec = dateVec, envType = 'woa')

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#
# create a list of all the desired input likelihood rasters
L.rasters1 <- list(L.sst = L.sst.par, L.light = L.light)
L.rasters2 <- list(L.sst = L.sst.par, L.light = L.light, L.ohc = L.ohc.par)
L.rasters3 <- list(L.sst = L.sst.par, L.light = L.light, L.prof = L.hycom.par)
L.rasters4 <- list(L.sst = L.sst.par, L.light = L.light, L.prof = L.qwoa.par)

# L.sst is the resolution/extent we're sampling everything TO
L.res1 <- resample.grid.par(L.rasters1, L.rasters1$L.sst)
L.res2 <- resample.grid.par(L.rasters2, L.rasters2$L.sst)
L.res3 <- resample.grid.par(L.rasters3, L.rasters3$L.sst)
L.res4 <- resample.grid.par(L.rasters4, L.rasters4$L.sst)

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
                    L2 = L.res2[[1]]$L.light,
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
par0 <- calc.param(migr.spd = 4, resid.spd = 2, g = g)
D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode
K1 <- gausskern(D1[1], D1[2], muadv = 0)
K2 <- gausskern(D2[1], D2[2], muadv = 0)

# RUN THE FILTER STEP
f <- hmm.filter(g, L, K1, K2, P.final1)

# plot if you want to see confidence limits
#res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
#fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P.final1)

# plot if you want to see confidence limits
#sres = apply(s[1,,,], 2:3, sum, na.rm=T)
#fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#

tr <- calc.track(s, g, dateVec)

plotHMM(s, tr, dateVec, ptt, save.plot = F)
#setwd(myDir)
write.table(tr, file=paste(ptt, '_HMM_track.csv', sep=''), sep = ',', col.names = T)
base::save.image(paste(ptt,'_hmm_geo.RData', sep=''))

#=======================================================================================#
## END
#=======================================================================================#



spot <- read.table('141259_crawl_track.csv',sep=',',header=T)
plot(tr$lon,tr$lat,type='l')
world(add=T,fill=T)
lines(spot$lon,spot$lat,col='blue')

pdf('check_141259_L.pdf',width=12,height=16)
par(mfrow=c(2,2))
for (i in 1:length(dateVec)){
  plot(L.res2[[1]]$L.light[[i]]); world(add=T); points(spot$lon[i], spot$lat[i]); title(dateVec[i])
  plot(L.res2[[1]]$L.sst[[i]]); world(add=T); points(spot$lon[i], spot$lat[i])
  plot(L.res2[[1]]$L.ohc[[i]]); world(add=T); points(spot$lon[i], spot$lat[i])
  image.plot(lon, lat, L[i,,]); world(add=T); points(spot$lon[i], spot$lat[i], pch=16, col='white'); title(dateVec[i])
  #image.plot(lon, lat, f$pred[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  #image.plot(lon, lat, s[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  #image.plot(lon, lat, s[2,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
}
dev.off()

pdf('check_141259_f_s.pdf',width=12,height=16)
par(mfrow=c(2,2))
for (i in 1:length(dateVec)){
  image.plot(lon, lat, L[i,,]); world(add=T); points(spot$lon[i], spot$lat[i], col='white'); title(dateVec[i])
  image.plot(lon, lat, f$pred[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  image.plot(lon, lat, s[1,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
  image.plot(lon, lat, s[2,i,,]); world(add=T); points(spot$lon[i], spot$lat[i])
}
dev.off()

