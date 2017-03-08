#=========================
#modify example from ws2trackit.r to fit one of our sharks
#=========================
#ensure LightLoc.csv file has begin data in row 4 and end data in row 6 with inittime (format as custom mm/dd/yy hh:mm), initlat, and initlon filled in.
#load wc2trackit & fix.trackit.dates functions from C:/...Documents/RCode/wc2trackit.r  **NOT** wc2trackit_dev.r
library(trackit)
library(kftrack)
library(ukfsst)
wd <- 'C:/RData/HMMoce_run/data/141259/'
setwd(wd)

lfile <- '141259-LightLoc.csv'
track <- wc2trackit(lfile, csv=T)
sp = track$sp; sp[3] <- sp[3]+2000
ep = track$ep; ep[3] <- ep[3]+2000
#track = fix.trackit.dates(track)
#track = two.layer.depth.corr(track, daybyday = FALSE, D0 = 50)
track <- as.data.frame(track$tdata)
names(track) <- list('year','month','day','hour','min','sec','depth','light','temp')
ptrack = prepit(track, fix.first = sp, fix.last=ep, scan=F)
fit = trackit(ptrack)
fit2 = trackit(ptrack, D.ph=2)
fit3 = trackit(ptrack, D.ph=-1, D.init=300)
fit4= trackit(ptrack, D.ph=-1, D.init=200)
#fit6=trackit(ptrack, ss3.ph=4)
#fit7=trackit(ptrack, ss3.ph=-1)
#fit9=trackit(ptrack, ss1.ph=-1)
#fit11= trackit(ptrack, D.ph=-1, D.init=300, ss3.ph=-1, ss2.ph=-1, ss1.ph=-1)

fitmap(fit4)

setwd(wd)
combine = data.frame(cbind(fit4$most.prob.track, fit4$date))
combine[,3] = combine[,3] + 60 * 365.25
write.table(combine, file="141259_trackit.csv", sep=",", quote=FALSE, append=FALSE, na="NA")

#===================
# Bathymetric Correction
#===================

setwd(wd)
bath = get.bath.data(-75, -30, 10, 55, folder = tempdir(), seaonly = TRUE, res = 1)

bath$data[bath$data > -15] = 1

setwd(wd)
dts <- as.POSIXct(fit4$date, format='%d%b%Y')
xtrack <- as.data.frame(cbind(day(dts), month(dts), year(dts), fit4$most.prob.track[,1:2]))
xtrack[,6:7] <- NA
names(xtrack) <- list('Day','Month','Year','Lon','Lat','SST','Max_depth')

maxd <- read.table('141259-MinMaxDepth.csv', sep=',', header=T, blank.lines.skip=F)
maxd <- maxd[,c(5,12)]
depdts <- as.POSIXct(maxd$Date, format='%H:%M:%S %d-%b-%Y')
depdts <- as.Date(depdts)

for (i in 1:length(xtrack[,1])){
  xtrack[i,7] <- maxd[which(depdts %in% as.Date(dts[i])),2]
}

fmat = prepb(fit4, xtrack)

fmat$Lon_E = fmat$Lon_E - 360 #reformats fmat so bathy correct actually works with correct long and depth
fmat$max_depth = fmat$max_depth * -1
str(fmat)

btrack = make.btrack(fmat,bath)

write.table(btrack, file="141259_fit_btrack.csv", sep=",", quote=FALSE, append=FALSE, na="NA", col.names=TRUE)

#btrack_fitmap=as.data.frame(matrix(0,length(btrack[,1]),2))
#btrack_fitmap[,1]=btrack$Lon_E
#btrack_fitmap[,2]=btrack$Lat_N
#btrack_fitmap

#plotbasemap(-75,-30,10,55)
#image(bath[[1]],bath[[2]],t(bath[[3]]),col=bath.colors(100), zlim=c(-5000,0), xlim=c(36,43), ylim=c(15,22), xlab="Longitude", ylab="Latitude", add=T)
#lines(btrack_fitmap, col="blue")
#points(btrack_fitmap, pch=19)

save.image("141259_traditional_geo.RData")


#========================================================#
#export locations from GPE derived points and run them through kftrack

setwd(wd)
lfile='141259-Locations-GPE2.csv'

data = read.table(lfile,sep=",",header=TRUE)
data = data[,c(2,4,5:14)]
data.gpe = data[data[,3]=='GPE',]
plotbasemap(min(data.gpe$Longitude)-5,max(data.gpe$Longitude)+5,min(data.gpe$Latitude)-5,max(data.gpe$Latitude)+5)
points(data.gpe$Longitude,data.gpe$Latitude)

dformat <- '%Y-%m-%d %H:%M:%S'
ddates = as.POSIXct(strptime(as.character(data.gpe$Date),format = dformat)) #reads dates as dates
ddates = as.POSIXct(trunc(ddates, 'day'))
year = as.numeric(format(ddates, '%Y')) #extracts year
month = as.numeric(format(ddates, '%m')) #extracts month
day = as.numeric(format(ddates, '%d')) #extracts day of month

blankline = c(0,0,0,0,0)
data.kf = cbind(day, month, year, data.gpe[,c(6,5)])
data.kf = rbind(blankline, data.kf, blankline)
data.kf[1,] = c(sp[c(5,4,3,1,2)])
data.kf[length(data.kf[,1]),] = c(ep[c(5,4,3,1,2)])
#data.kf = data.kf[data.kf$Latitude >= 5,]

kfit = kftrack(data.kf)
#kfit1 = kftrack(data.kf, D.a=F, D.i=200)
#kfit2 = kftrack(data.kf,D.a=F,D.i=500,bx.a=F,by.a=F,sx.a=F)
kfit3=kftrack(data.kf,D.a=F,bx.a=F,by.a=F,sy.a=F,sx.a=F)
kfit4=kftrack(data.kf,D.a=F,bx.a=F,sx.a=F)
kfit5=kftrack(data.kf,D.a=F,D.i=300,bx.a=F,by.a=F,sx.a=F)
kfit6=kftrack(data.kf,D.a=F,D.i=300,bx.a=F,by.a=F,sy.a=F,sx.a=F)
kfit7=kftrack(data.kf,D.a=F,D.i=300,bx.a=F,sx.a=F)
kfit8=kftrack(data.kf,D.a=F,D.i=200,bx.a=F,sx.a=F)
kfit9=kftrack(data.kf,D.a=F,D.i=200,v.a=F)
#kfit10=kftrack(data.kf,D.a=F,D.i=200,u.a=F,v.a=F)
kfit11=kftrack(data.kf,D.a=F,D.i=200,u.a=F,v.a=F,bx.a=F)
kfit12=kftrack(data.kf,D.a=F,D.i=200,u.a=F,v.a=F,bx.a=F,by.a=F)
kfit13=kftrack(data.kf,D.a=F,u.a=F,v.a=F,bx.a=F)
kfit14=kftrack(data.kf,u.a=F,v.a=F,bx.a=F)
kfit = kfit11

plotbasemap(min(kfit$most.prob.track[,1])-5,max(kfit$most.prob.track[,1])+5,min(kfit$most.prob.track[,2])-5,max(kfit$most.prob.track[,2])+5)
#points(kfit$most.prob.track[,1],kfit$most.prob.track[,2])
lines(kfit$most.prob.track[,1],kfit$most.prob.track[,2])

#write.table(kfit$most.prob.track, file="141259_kfit.csv", sep=",", quote=FALSE, append=FALSE, na="NA")
#write.table(kfit$date, file="141259_kfitdate.csv", sep=",", quote=FALSE, append=FALSE, na="NA")

setwd(wd)
maxdepth=read.table(file='141259-MinMaxDepth.csv',sep=',',blank.lines.skip=FALSE,header=TRUE)
maxdepth=maxdepth[,c(2,5,12)]
dates = as.POSIXct(strptime(as.character(maxdepth$Date),format = '%H:%M:%S %d-%b-%Y')) #reads dates as dates
#times = as.POSIXct(as.character(maxdepth$Time),format = '%H:%M:%S')
month=as.numeric(format(dates,'%m'))
day=as.numeric(format(dates,'%d'))
year=as.numeric(format(dates,'%Y'))
maxdepth=cbind(maxdepth,month,day,year)
maxdepth[,7]=paste(maxdepth$month,'/',maxdepth$day,'/',maxdepth$year,sep='')
colnames(maxdepth)=list("ptt","datetime","maxdepth","month","day","year","date")
kfdates=paste(data.kf[,2],'/',data.kf[,1],'/',data.kf[,3],sep='')

for (i in 1:length(kfdates)){
  data.kf[i,7]=max(maxdepth[which(maxdepth[,7]==kfdates[i]),3],na.rm=TRUE)
}

str(data.kf)
colnames(data.kf)=list('day','month','year','long','lat','sst','maxdepth')
# write.table(data.kf,file='141259_data.csv',sep=',',row.names=FALSE,col.names=TRUE)

data.kf$maxdepth=data.kf$maxdepth*-1
fmat = prepb(kfit, data.kf)
str(fmat)
summary(fmat)
#fmat$Lon_E=fmat$Lon_E+360 #make sure longs 0-360

bath = get.bath.data(min(fmat$Lon_E)-5,max(fmat$Lon_E)+5,min(fmat$Lat_N)-5,max(fmat$Lat_N)+5, folder = tempdir(), seaonly = TRUE, res = .5)
bath$data[bath$data>-15]=1

kfbtrack = make.btrack(fmat,bath)

btrack = kfbtrack
#btrack[,c(12:13)]=kfit11$SST[,2:3]
colnames(btrack)=list('year','month','day','v11','v12','v21','v22','long','lat','maxdepth','o_sst')
getwd()
write.table(btrack, file="141259_kfit11_btrack.csv", sep=",", quote=FALSE, append=FALSE, na="NA", col.names=TRUE)

fitmap(kfit)
lines(btrack[,8:9])
save.image('141259_traditional_geo.RData')
