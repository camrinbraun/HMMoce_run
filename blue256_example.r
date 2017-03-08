# RUN BLUE SHARK EXAMPLE, 141256, USING HMMoce
library(HMMoce)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141256/') 

# READ IN TAG DATA
ptt <- '141256'

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                              24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
myDir <- '~/Documents/WHOI/RCode/HMMoce/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
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
locs <- read.table('141259-Locations-GPE2.csv', sep = ',', header = T, blank.lines.skip = F)

#----------------------------------------------------------------------------------#
# FURTHER PREPARATION
# Set spatial limits and download env data
#----------------------------------------------------------------------------------#

# SET SPATIAL LIMITS, IF DESIRED
sp.lim <- list(lonmin = -95, lonmax = -52, latmin = 10, latmax = 55)

if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}

# IF USING SST, DOWNLOAD THE SST DATA BY UNCOMMENTING THE GET.ENV FUNCTION
sst.dir <- paste('~/Documents/WHOI/RData/SST/OI/', ptt, '/',sep = '')
#get.env(sst.udates[1], type = 'sst', spatLim = sp.lim, save.dir = sst.dir)

# IF USING OHC, DOWNLOAD HYCOM DATA BY UNCOMMENTING THE GET.ENV FUNCTION
hycom.dir <- paste('~/Documents/WHOI/RData/HYCOM/', ptt, '/',sep = '')
#get.env(pdt.udates, type = 'ohc', spatLim = sp.lim, save.dir = hycom.dir)

#----------------------------------------------------------------------------------#
# CALC LIKELIHOODS
#----------------------------------------------------------------------------------#

# LIGHT LIKELIHOOD
#L.light <- calc.light(light, locs.grid = locs.grid, dateVec = dateVec)
# OR
L.light <- calc.locs(locs, iniloc = iniloc, locs.grid = locs.grid, dateVec = dateVec, errEll = TRUE, gpeOnly = TRUE)
L.light <- L.light$L.locs

#-------
# GENERATE DAILY SST LIKELIHOODS
L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec)

#-------
# GENERATE DAILY OCEAN HEAT CONTENT (OHC) LIKELIHOODS
L.ohc <- calc.ohc(pdt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '')

#-------
# GENERATE DAILY PROFILE LIKELIHOODS
L.prof <- calc.profile(pdt, dateVec = dateVec, envType = 'hycom', hycom.dir = hycom.dir)
#L.prof.woa <- calc.profile(pdt, dat = woa, lat = lat, lon = lon, dateVec = dateVec, envType = 'woa')

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#

L.rasters <- list(L.ohc = L.ohc, L.sst = L.sst, L.light = L.light)
L.res <- resample.grid(L.rasters, L.rasters$L.sst)

L.mle.res <- L.res$L.mle.res
g <- L.res$g; lon <- g$lon[1,]; lat <- g$lat[,1]
g.mle <- L.res$g.mle

#----------------------------------------------------------------------------------#
# LOAD AND FORMAT DATAFRAME OF KNOWN LOCATIONS, IF ANY
#----------------------------------------------------------------------------------#

#colnames(known.locs) <- list('date','lat','lon')
#   where 'date' is from as.Date(known.locs$date)

#----------------------------------------------------------------------------------#
# COMBINE LIKELIHOOD MATRICES
#----------------------------------------------------------------------------------#

L <- make.L(L1 = L.res[[1]]$L.ohc, L2 = L.res[[1]]$L.sst, 
            L3 = L.res[[1]]$L.light,
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc)

L.mle <- L$L.mle; L <- L$L

#----------------------------------------------------------------------------------#
# TRY THE MLE.

# NOT RIGHT NOW
#t1 <- Sys.time()
#guess <- c(log(10),log(10),log(0.5),log(0.5),log(0.95/0.05),log(0.95/0.05))
#fit <- nlm(neg.log.lik.fun, guess, g.mle, L.mle)
#t2 <- Sys.time()

## **THESE OUTPUT PARAMETERS ARE PIXEL-BASED. DON'T FORGET TO CONVERT FOR USE
##  WITH THE HIGHER RESOLUTION LIKELIHOOD RESULTS STORED IN L 
D1 <- exp(fit$estimate[1:2]) # parameters for kernel 1. this is behavior mode transit. log-transformed movement parameters (diffusivities) pertaining
# to the first behavioural state.

D2 <- exp(fit$estimate[3:4]) # parameters for kernel 2. resident behavior mode. log-transformed movement 
# parameters (diffusivities) pertaining to the second behavioural state.

p <- 1/(1+exp(-fit$estimate[5:6])) # logit-transformed
#transition probabilities for switching between the two behavioural states 
#Probably need to express kernel movement in terms of pixels per time step.
#The sparse matrix work likely renders this unnecessary, but going back to 
#gausskern, it is. For example, if we have .25 degree and daily time step,
#what would the speed of the fish be when moving fast? 4 pixels/day?

#----------------------------------------------------------------------------------#
# OR... JUST DEFINE THE PARAMETERS
par0=c(8.908,10.27,1.152,0.0472,0.707,0.866)
D1 <- par0[1:2] # parameters for kernel 1. this is behavior mode transit
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode
p <- par0[5:6]

#----------------------------------------------------------------------------------#
# GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
K1 = gausskern(D1[1], D1[2], muadv = 0)
K2 = gausskern(D2[1], D2[2], muadv = 0)
P <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)

#----------------------------------------------------------------------------------#
# RUN THE FILTER STEP
f = hmm.filter(g, L, K1, K2, P)

# plot if you want to see confidence limits
#res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
#fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P)

# plot if you want to see confidence limits
#sres = apply(s[1,,,], 2:3, sum, na.rm=T)
#fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#
T <- dim(s)[2]
meanlat <- apply(apply(s, c(2, 4), sum) * repmat(t(as.matrix(g.mle$lat[,1])), T, 1), 1, sum)
meanlon <- apply(apply(s, c(2, 3), sum) * repmat(t(as.matrix(g.mle$lon[1,])), T, 1), 1, sum)

#**track <- calc.track(distr, g)**

plot(meanlon,meanlat,type='l')
world(add=T, fill=T, col='grey')

#=======================================================================================#
## END
#=======================================================================================#
