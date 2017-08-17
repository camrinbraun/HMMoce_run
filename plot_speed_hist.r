
require(trip)
#locs <- read.table('~/Documents/WHOI/Data/Blues/2015/141254/141270-Locations.csv', sep=',', header=T)
#locs <- read.table('~/Documents/WHOI/Data/Blues/2015/141256/141268-Locations.csv', sep=',', header=T)
#locs <- read.table('~/Documents/WHOI/Data/Blues/2015/141259/141261-Locations.csv', sep=',', header=T)
locs <- read.table('~/Documents/WHOI/Data/Makos/2015/141257/141267-Locations.csv', sep=',', header=T)

nameList <- list('Ptt','Date','Quality','Latitude','Longitude')
locs <- locs[,which(colnames(locs) %in% nameList)]

# get rid of Z pings
locs <- locs[which(locs$Quality != 'Z'),]


if(class(locs$Date) == 'factor'){
  locs$Date <- as.POSIXct(locs$Date, format = HMMoce::findDateFormat(locs$Date), tz = 'UTC')
}

locs <- locs[which(!duplicated(locs$Date)),]

track <- locs

track <- track[order(track$Date),]

# convert track to SpatialPointsDataFrame
coordinates(track) <- ~Longitude + Latitude
proj4string(track)=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# then to trip object
tr <- trip(track, c('Date','Ptt'))

tdist <- trip::trackDistance(tr, longlat = T)
timediff <- diff(tr@data$Date) # in seconds

if(units(timediff) == 'hours') timediff <- timediff * 3600
if(units(timediff) == 'mins') timediff <- timediff * 60
#if(units(timediff) == 'seconds') timediff <- unclass(timediff)

timediff <- unclass(timediff)

spd <- (tdist[2:length(tdist)] / timediff) * 1000 # m/s

#h <- hist(spd, density=T, xlab = 'Speed (m/s)', main='', xlim=c(0,10), breaks=seq(0,ceiling(max(spd)*1000), by=.25))#, ylim=c(0,50))

d257 <- density(spd, na.rm=T, from=0, to=10)

plot(d254, xlim=c(0,4), ylim=c(0,1), main='', xlab='Speed (m/s)')
lines(d256, col='red')
lines(d257, col='green')
lines(d259, col='blue')

df254 <- data.frame(cbind(d254$x, d254$y))
names(df254) <- c('x','y')
df254$group <- 254

df256 <- data.frame(cbind(d256$x, d256$y))
names(df256) <- c('x','y')
df256$group <- 256

df257 <- data.frame(cbind(d257$x, d257$y))
names(df257) <- c('x','y')
df257$group <- 257

df259 <- data.frame(cbind(d259$x, d259$y))
names(df259) <- c('x','y')
df259$group <- 259

df <- rbind(df254, df256, df257, df259)
df$group[which(df$group == 254)] <- 141254
df$group[which(df$group == 256)] <- 141256
df$group[which(df$group == 257)] <- 141257
df$group[which(df$group == 259)] <- 141259

ggplot(df, aes(x=x, y=y, col=as.factor(group))) + geom_line(lwd=.75) + xlim(0,4) +
  xlab('Speed (m/s)') + ylab('Density') + 
  scale_colour_manual(values=c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3'), name='PTT')

#==========================

hist(spd*1000, density=T, xlab = 'Speed (m/s)', main='', xlim=c(0,10), breaks=seq(0,ceiling(max(spd)*1000), by=.25))#, ylim=c(0,50))
#hist(spd*1000, xlab = 'Speed (m/s)', main='')#, xlim=c(0,30), breaks=seq(0,ceiling(max(spd)), by=1), ylim=c(0,50))
#abline(v=mean(h[[3]]), col='red')


# Compare MPG distributions for cars with 
# 4,6, or 8 cylinders
library(sm)

# create value labels 
cyl.f <- factor(cyl, levels= c(4,6,8),
                labels = c("4 cylinder", "6 cylinder", "8 cylinder")) 

# plot densities 
sm.density.compare(mpg, cyl, xlab="Miles Per Gallon")
title(main="MPG Distribution by Car Cylinders")


