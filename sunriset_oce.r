# test sunrise optimization with MTI SRSS data

library(rgdal)
library(oce)
library(RODBC)
library(maptools)
library(maps)
library(lubridate)

# start point
sp = c(-64.700, 43.540)

# OCE methods
# rise <- as.POSIXct("2014-10-17 10:19:00", tz="UTC") #+ 4*3600
# set <- as.POSIXct("2014-10-17 21:44:00", tz="UTC") #+ 4*3600
spm = t(matrix(sp))
sday = as.POSIXct("2014-6-17", tz = 'UTC')

rise <- sunriset(spm,  sday, direction = "sunrise", POSIXct.out=TRUE)$time
set <- sunriset(spm, sday, direction = 'sunset', POSIXct.out = TRUE)$time

mismatch <- function(lonlat) 
{
  sunAngle(rise, lonlat[1], lonlat[2])$altitude^2 + sunAngle(set, lonlat[1], lonlat[2])$altitude^2
}
result <- optim(c(1,1), mismatch)

dist <- geodDist(result$par[1], result$par[2], sp[1], sp[2])
cat(sprintf("Infer start point latitude %.2f and longitude %.2f; distance mismatch %.0f km", 
            result$par[2], result$par[1], dist))

msp = as.matrix(t(sp), nrow=1, ncol = 2)
ssp = SpatialPoints(msp, proj4string=CRS("+proj=longlat +datum=WGS84"))

sdate <- as.POSIXct("2014-10-17", tz="UTC") #+ 4*3600

# correct sunrise time in UTC
crise = sunriset(msp, sdate, direction = 'sunrise', POSIXct.out=TRUE)[2]

# correct sunset time in UTC
cset = sunriset(msp, sdate, direction = 'sunset', POSIXct.out=TRUE)[2]

# how much off is tre first day at liberty from the actual start position
crise[1,1] - rise
# Time difference of 16.27241 mins

cset[1,1] - set
# Time difference of -11.69904 mins

# maptools methods
# look at a matrix of positions
lonvec = seq(sp[1]-5, sp[1]+5, by=.25)
latvec = seq(sp[2]-5, sp[2]+5, by=.25)

llmat = expand.grid(lonvec, latvec)

llmat = SpatialPoints(llmat, proj4string=CRS("+proj=longlat +datum=WGS84"))

srmat = sunriset(llmat, sdate, direction = 'sunrise', POSIXct.out=TRUE)[2]
ssmat = sunriset(llmat, sdate, direction = 'sunset', POSIXct.out=TRUE)[2]

# time difference for sr and ss by grid cell
srdiff = srmat[,1]-rise
ssdiff = ssmat[,1]-set

llmat$srdiff = abs(as.numeric(srdiff))
llmat$ssdiff = abs(as.numeric(ssdiff))


# plot differences 
par(mfrow=c(1,2))
# plot(srmat, pch = 19, cex.=.03, col = 4)
# plot(ssmat, pch = 19, cex.=.03, col=3)


plot(llmat, cex = llmat$srdiff/max(llmat$srdiff), pch=1)
plot(ssp, add=T, col=2, pch=19, cex=1.4)
degAxis(1)
degAxis(2)
title('sunrise difference')
legend('bottom', , legend = round(quantile(llmat$srdiff/60)[2:5]), horiz = T, pch = 1, pt.cex = quantile(llmat$srdiff/max(llmat$srdiff))[2:5], title = 'minutes')
map('world', add=T)

plot(llmat, cex = llmat$ssdiff/max(llmat$ssdiff), pch=1)
plot(ssp, add=T, col=2, pch=19, cex=1.4)
degAxis(1)
degAxis(2)
title('sunset difference')
legend('bottom', , legend = round(quantile(llmat$srdiff/60)[2:5]), horiz = T, pch = 1, pt.cex = quantile(llmat$srdiff/max(llmat$srdiff))[2:5], title = 'minutes')
map('world', add=T)


# look at one point for a year

dates = seq(ymd('2010-1-1', tz= 'UTC'), ymd('2010-12-31', tz= 'UTC'), by = 'day')
newlon = numeric(length = 365)
for(i in 1:365){
  
  spm = t(matrix(sp))
  rise <- sunriset(spm,  dates[i], direction = "sunrise", POSIXct.out=TRUE)$time
  set <- sunriset(spm, dates[i], direction = 'sunset', POSIXct.out = TRUE)$time
  
  mismatch <- function(lonlat) 
  {
    sunAngle(rise, lonlat[1], lonlat[2])$altitude^2 + sunAngle(set, lonlat[1], lonlat[2])$altitude^2
  }
  
  newlon[i] <- optim(c(1,1), mismatch)$par[1]
  
}

plot(dates, newlon)
abline(h = sp[1], col = 2, lty = 2, lwd = 2)
title('one year of optimized longitudes \n with known sr and ss')


# now, train the longitude based on tag measured sunrise and sunset and the pervious days position... 

ll = read.csv('C:/Users/ben/Downloads/141256-LightLoc.csv', skip = 2, header = T)

ll = ll[4:nrow(ll),]

sp = t(as.matrix(c(-69.42263, 41.5750)))
ep = t(matrix(c(-69.0986, 26.6166)))


# add the CRAWL track

ctrack = read.csv('C:/Users/ben/Downloads/141256_crawl_track.csv')
ctrack$date = ymd(ctrack$date, tz ='UTC')

# add theGPE2 track

gtrack = read.csv('C:/Users/ben/Downloads/141256-Locations-GPE2.csv')
gtrack$Date = ymd_hms(gtrack$Date, tz ='UTC')



# add time from Delta to approximate midpoint of light measurement time

ll$dtime = dmy_hms(paste(ll$Day, ll$Time), tz = 'UTC')+ll$Delta*4

ll$date = as.POSIXct(trunc(ll$dtime, 'day'))

dawnidx  = ll$Type=='Dawn'
duskidx  = ll$Type=='Dusk'

# add and subtract time to sunset and sunrise to account for depth
# one minute every 10 meters??
# use Kd??
IO = function(IZ, K, z) IZ/exp(-K*z)

ll$middepth = (ll$MaxDepth-ll$MinDepth)/2


ll$mean_light = rowMeans(ll[,grep('LL', names(ll))[1:8]])

# sunrise and sunset times

rise.v = ll$dtime[dawnidx]
set.v = ll$dtime[duskidx]

# adjustment
ll$adj = 60*2.5*log(abs(ll$mean_light-IO(ll$mean_light, .1, ll$middepth))) # adjustment in minutes
ll$adj[is.infinite(ll$adj)] = 0


rise.v = rise.v-ll$adj[dawnidx]
set.v = set.v +ll$adj[duskidx]

newlon = newlon2 = newlat = numeric(length = length(rise.v))
newlon[1] = newlon2[1] = sp[1]

for(i in 2:length(rise.v)){

  rise <- rise.v[i]
  set <- set.v[i]
  
  mismatch <- function(sp) 
  {
    sunAngle(rise, sp[1], sp[2])$altitude^2 + sunAngle(set, sp[1], sp[2])$altitude^2
  }
  # 
  mismatch2 <- function(sp)
  {
    sunAngle(rise, sp[1], sp[2])$altitude^2 + sunAngle(set, sp[1], sp[2])$altitude^2 + sunAngle(rise, ep[1], ep[2])$altitude^2 + sunAngle(set, ep[1], ep[2])$altitude^2
  }
  
  if (i >=5) {
    medlon = median(newlon[(i-4):i])
    lim = 2*sd(newlon[(i-4):i])
  }else{
  medlon = newlon[i-1]  
  lim = 2  
  }
  
  
  res = optim(c(sp[1] ,sp[2]), mismatch, lower = c(medlon-lim, sp[2]-5), upper = c(medlon+lim, sp[2]+5) )
  
  # res = optim(c(sp[1] ,sp[2]), mismatch)
  
  res2 =   res = optim(c(ep[1] ,ep[2]), mismatch, lower = c(newlon2[i-1]-2, ep[2]-5), upper = c(newlon2[i-1]+2, ep[2]+5) )
  # res2 =   res = optim(c(ep[1] ,ep[2]), mismatch, lower = c(newlon[i-1]-3, ep[2]-10), upper = c(newlon[i-1]+3, ep[2]+10) )
  
  # res = optim(c(sp[1] ,sp[2]), mismatch, lower = c(newlon[i-1]-3, sp[2]-10), upper = c(newlon[i-1]+3, sp[2]+10) )
  
  # res = optim(c(sp[1] ,sp[2]), mismatch)
  
  # res2 =   res = optim(c(ep[1] ,ep[2]), mismatch2, lower = c(newlon2[i-1]-3, ep[2]-10), upper = c(newlon2[i-1]+3, ep[2]+10) )
  
  newlon[i] <- res$par[1]
  newlon2[i] = res2$par[1]
}

##
par(mfrow=c(2,1))

idx = newlon > quantile(newlon, .75)|newlon < quantile(newlon, .25)

plot(rise.v, newlon, pch = 19)
points(rise.v[idx], newlon[idx], col = 2, pch = 19)
title('flag middle 50% quantile')

# plot(rise.v[!idx], newlon[!idx], col = 1, pch = 19, lty = 2, typ = 'o')
plot(rise.v, newlon, col = 1, pch = 19, lty = 2, typ = 'o')
lines(ctrack$date, ctrack$lon, col = 4, lwd = 2)
lines(gtrack$Date, gtrack$Longitude, col = 3, lwd = 2)
grid()

lines(rise.v, newlon2, col = 2)


# par(mfrow=c(2,1))
# plot(rise.v, newlon, typ='l')
# plot(rise.v, ll$MaxDepth[ll$Type=='Dawn']*-1, typ='l', col = 2, axes = T)

# 
# par(mfrow=c(2,1))
# 
# plot(rise.v, newlon, typ='l')
# lines(ctrack$date, ctrack$lon, col = 4, lwd = 2)
# grid()
# 
# plot(rise.v, ll$MaxDepth[ll$Type=='Dawn']*-1, typ='l', col = 2, axes = T, ylab = 'Max Depth')
# grid()
# 
# # plot quantiles of newlon
# hist(newlon, breaks = 25)
# 
# abline(v = quantile(newlon), col = 2, lwd = 2, lty = 2)
# 
# idx = newlon > quantile(newlon, .75)|newlon < quantile(newlon, .25)

# 
# # now iterate each day
# 
# newlon = numeric(length = length(rise.v))
# for(i in 1:length(rise.v)){
#   
#   spm = t(matrix(sp))
#   rise <- rise.v[i]
#   set <- set.v[i]
#   
#   mismatch <- function(lonlat) 
#   {
#     sunAngle(rise, lonlat[1], lonlat[2])$altitude^2 + sunAngle(set, lonlat[1], lonlat[2])$altitude^2
#   }
#   
#   if(i == 1 ){
#     res = optim(c(sp[1] ,sp[2]), mismatch)
#     newlon[i] <- res$par[1]
#   }else{
#      res = optim(c(res$par[1], res$par[2]), mismatch)
#      newlon[i] <- res$par[1]
#   }
# }
# 
# 

#--------------------------------------------------------------------#
librry(plyr)
ddepth = ddply(ll, 'date', function(x) x$MaxDepth[1] - x$MinDepth[1])
ddepth$newlon = newlon[2:(length(newlon)-1)]

ddepth = merge(ddepth, ctrack, by = 'date')

plot(ddepth$V1, ddepth$lon - ddepth$newlon)

