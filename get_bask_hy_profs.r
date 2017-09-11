dataDir <- '~/ebs/Data/BaskingSharks/batch/'
library(raster); library(rgeos); library(RNetCDF)
meta <- read.table(paste(dataDir, 'bask_metadata_v2.csv',sep=''), sep=',', header=T)
meta$PopDate <- as.Date(meta$PopDate, format='%m/%d/%y')
meta <- meta[which(!is.na(meta$res)),]

profs <- list()
for (tt in 1:nrow(meta)){
  ptt <- meta$PTT[tt]
  dte <- meta$PopDate[tt]
  if (meta$bnds[tt] == 'big'){
    ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/big/bask_big_', dte, '.nc', sep='')
  } else{
    ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/small/bask_', dte, '.nc', sep='')
  }
  
  if(!file.exists(ncName)){
    if (meta$bnds[tt] == 'big'){
      ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/big/bask_big_', dte-1, '.nc', sep='')
    } else{
      ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/small/bask_', dte-1, '.nc', sep='')
    }
  }
  
  if(!file.exists(ncName)){
    sp.lim <- list(lonmin=meta$PopLong[tt]-1, lonmax=meta$PopLong[tt]+1,
                   latmin=meta$PopLat[tt]-1, latmax=meta$PopLat[tt]+1)
    if (meta$bnds[tt] == 'big'){
      hycom.dir <- '~/ebs/EnvData/hycom3/BaskingSharks/big/'
      get.env(dte, filename='bask_big', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
      ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/big/bask_big_', dte, '.nc', sep='')
      
    } else{
      hycom.dir <- '~/ebs/EnvData/hycom3/BaskingSharks/small/'
      get.env(dte, filename='bask', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
      ncName <- paste('~/ebs/EnvData/hycom3/BaskingSharks/small/bask_', dte, '.nc', sep='')
      
    }
    
  }
  
  # open nc and get the indices for the vars
  nc1 =  RNetCDF::open.nc(ncName)
  ncnames = NULL
  nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
  for(ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1, ii)$name
  temp.idx <- grep('temp', ncnames, ignore.case=TRUE) - 1
  lat.idx <- grep('lat', ncnames, ignore.case=TRUE) - 1
  lon.idx <- grep('lon', ncnames, ignore.case=TRUE) - 1
  dep.idx <- grep('dep', ncnames, ignore.case=TRUE) - 1
  
  # get attributes, if they exist
  ncatts <- NULL
  nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
  for(ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1, temp.idx, ii)$name
  scale.idx <- grep('scale', ncatts, ignore.case=TRUE) - 1
  if(length(scale.idx) != 0){
    scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=scale.idx)
  } else{
    scale <- 1
  }
  off.idx <- grep('off', ncatts, ignore.case=TRUE) - 1
  if(length(off.idx) != 0){
    offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute=off.idx)
  } else{
    offset <- 1
  }
  
  # get and check the vars
  depth <- RNetCDF::var.get.nc(nc1, dep.idx)
  lon <- RNetCDF::var.get.nc(nc1, lon.idx)
  if(length(dim(lon)) == 2) lon <- lon[,1]
  if(!any(lon < 180)) lon <- lon - 360
  lat <- RNetCDF::var.get.nc(nc1, lat.idx)
  if(length(dim(lat)) == 2) lat <- lat[1,]
  
  # position indices
  elo <- which.min(abs(lon - meta$PopLong[tt]))
  ela <- which.min(abs(lat - meta$PopLat[tt]))
  
  # get profile
  temp <- var.get.nc(nc1, 'water_temp', start=c(elo,ela,1,1), count=c(1,1,NA,NA)) * scale + offset
  temp.idx <- which(!is.na(temp))
  
  profs[[tt]] <- list(ptt=ptt, temp=temp[temp.idx], depth=depth[temp.idx], pos=c(meta$PopLong[tt], meta$PopLat[tt]))
}


for (ii in 1:nrow(meta)){
  
  dataDir <- '~/ebs/Data/BaskingSharks/batch/'
  ptt <- meta$PTT[ii]
  setwd(paste(dataDir, '/', ptt, sep=''))
  
  iniloc <- data.frame(matrix(c(meta$TagDay[ii], meta$TagMonth[ii], meta$TagYear[ii], meta$TagLat[ii], meta$TagLong[ii], 
                                meta$PopDay[ii], meta$PopMonth[ii], meta$PopYear[ii], meta$PopLat[ii], meta$PopLong[ii]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  
  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day')) 
  
  # READ IN DATA FROM WC FILES
  myDir <- paste(dataDir, ptt, '/', sep='')
  setwd(myDir)

  # depth-temp profile data
  pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
  pdt.udates <- pdt$udates; pdt <- pdt$data
  tail(pdt)
  print(meta$PopDate[ii])
  pdt$Date <- as.Date(pdt$Date, format = findDateFormat(pdt$Date))
  
  pdt.i <- pdt[which(pdt$Date %in% (as.Date(pop)-2)),c(2,5:7)]
  pdt.i
  
  profs[[ii]]$pdt.i <- pdt.i

}

# re-run these to get pdt.i
#run.idx <- which(unlist(lapply(profs, FUN=function(x) nrow(x$pdt.i) == 0)))
#  [1]  5  6  7  9 12 16 22 23 25 26 28 29 30 33 34 35 36 37

setwd('~/ebs/Data/BaskingSharks/batch/')
#save(profs, file='bask_profile_list.rda')
load('bask_profile_list.rda')

for (tt in 1:nrow(meta)){
  temp.idx <- which(!is.na(profs[[tt]]$temp) & profs[[tt]]$depth <= max(profs[[tt]]$pdt.i$Depth))
  
  plot(profs[[tt]]$temp[temp.idx], profs[[tt]]$depth[temp.idx], ylim=c(350,0), type='l')#c(max(profs[[tt]]$depth[temp.idx]),0), type='l')
  points(profs[[tt]]$pdt.i$MinTemp, profs[[tt]]$pdt.i$Depth)
  points(profs[[tt]]$pdt.i$MaxTemp, profs[[tt]]$pdt.i$Depth)
  
}

# example profiles betw hycom and tags
# use ii = 2,4,19,21,26,29,30,31,36

idx <- c(26,29,30,31)

minT <- min(unlist(lapply(profs, FUN=function(x) min(x$temp, na.rm=T))[idx]))
maxT <- max(unlist(lapply(profs, FUN=function(x) max(x$temp, na.rm=T))[idx]))
maxZ <- max(unlist(lapply(profs, FUN=function(x) max(x$pdt.i$Depth, na.rm=T))[idx])) * 1.2

poly.cols <- c(rgb(228/255,26/255,28/255,130/255), rgb(55/255,126/255,184/255,130/255),
          rgb(77/255,175/255,74/255,130/255), rgb(152/255,78/255,163/255,130/255))
line.cols <- c(rgb(228/255,26/255,28/255), rgb(55/255,126/255,184/255),
               rgb(77/255,175/255,74/255), rgb(152/255,78/255,163/255))
textList <- c('A','B','C','D')
profs.sel <- list()
for (b in 1:length(idx)){
  profs.sel[[b]] <- profs[[idx[b]]]
  profs.sel[[b]]$line.col <- line.cols[b]
  profs.sel[[b]]$poly.col <- poly.cols[b]
  profs.sel[[b]]$text <- textList[b]
}

# set lims
xlims=c(-83,-28); ylims=c(-12,47)
x.at <- pretty(xlims, 5)
x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
y.at <- pretty(ylims, 5)
y.labels <- c(parse(text=paste(abs(y.at[which(y.at < 0)]), "*degree~S", sep="")), parse(text=paste(y.at[which(y.at >= 0)], "*degree~N", sep="")))

pdf('fig_hycom_profiles.pdf', width=11, height=8)
par(mfrow=c(1,2))

plot(0,0, type='n', xlim=xlims, ylim=ylims, axes=F, xlab='',ylab='')
fields::world(add=T, fill=T, col='grey80', border='grey80')
points(meta$PopLong, meta$PopLat)
points(meta$PopLong[idx], meta$PopLat[idx], col=cols, pch=16, cex=1.5)
text(-69.5, 44, 'A', col=line.cols[1], font=2)
text(-30, -6, 'B', col=line.cols[2], font=2)
text(-71, 29.8, 'C', col=line.cols[3], font=2)
text(-55, 20, 'D', col=line.cols[4], font=2)
#text(meta$PopLong[idx], meta$PopLat[idx], labels=meta$PTT[idx], col=cols)
#text(meta$PopLong[-idx], meta$PopLat[-idx], labels=meta$PTT[-idx])
axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels)
box()


plot(0,0, type='n', ylim=c(1100,-40), xlim=c(2.5, maxT), xlab = expression("Temperature " ( degree*C)), ylab='Depth (m)')
#lapply(profs.sel, FUN=function(x) points(x$pdt.i$MinTemp, x$pdt.i$Depth, col=x$col, pch=16))
lapply(profs.sel, FUN=function(x) polygon(c(x$pdt.i$MinTemp, rev(x$pdt.i$MaxTemp)), c(x$pdt.i$Depth, rev(x$pdt.i$Depth)), col=x$poly.col, border=x$poly.col))
lapply(profs.sel, FUN=function(x) lines(x$temp, x$depth, col=x$line.col, lwd=2))
lapply(profs.sel, FUN=function(x) text(x$temp[1], -30, x$text, col=x$line.col, font=2))

dev.off()


# maybe: 1
# yes: 2,3,4,6,7,8
# no: 5
plot(meta$PopLong, meta$PopLat)
points(meta$PopLong[idx], meta$PopLat[idx],col=cols,pch=16)
text(meta$PopLong[idx], meta$PopLat[idx], labels=meta$PTT[idx], col=cols)
#text(meta$PopLong[-idx], meta$PopLat[-idx], labels=meta$PTT[-idx])
world(add=T)

32 or 6
110493 OR 52557
9 or 5
52562 OR 52556

save.image('hycom_profile_fig_v1.rda')


idx <- c(19, 29, 31, 32)

32 or 6
