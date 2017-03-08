light <- read.table('141256-LightLoc.csv', sep=',', header=T, skip=2)

light.m <- light[,c('Day','Time','Type')]
light.m[,4] <- as.POSIXct(paste(light.m$Day, light.m$Time), format='%d-%b-%y %H:%M:%S', tz='UTC')
light.m[,5] <- NA
light.m[c(1:(nrow(light.m)-1)),5] <- paste(light.m[2:nrow(light.m),4])
light.m <- light.m[c(4:(nrow(light.m)-1)),]
light.m <- light.m[,c(4,5,3)]
light.m[,2] <- as.POSIXct(light.m[,2], tz='UTC')
light.m[which(light.m[,3] == 'Dawn'),4] <- 1
light.m[which(light.m[,3] == 'Dusk'),4] <- 2
light.m <- light.m[,c(1,2,4)]
names(light.m) <- list('tFirst', 'tSecond', 'type')

twl <- light.m
crds0 <- coord(twl, tol = 0)



#----------------
data(hoopoe2)
hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
residency <- with(hoopoe2, changeLight(tFirst,tSecond,type, rise.prob=0.1,
                                       set.prob=0.1, plot=FALSE, summary=FALSE))
hec <- HillEkstromCalib(hoopoe2,site = residency$site)



