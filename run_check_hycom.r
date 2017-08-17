# input spot data, hycom directory, pdt, ohc L, hycom L
require(HMMoce)
setwd('~/HMMoce_run/data/141254/')
load('141254_likelihoods_knock.RData')
L.ohc <- L.3; L.hycom <- L.5; L.woa <- L.4
hycom.dir <- '~/HMMoce_run/env_data/141254/hycom/'
spot <- read.table('141254_crawl_track.csv', sep=',', header=T)
str(pdt)

xlims <- c(-85, -46); ylims <- c(31, 47)
ex <- raster::extent(c(xlims, ylims))
L.ohc <- raster::crop(L.ohc, ex)

pdt$Date <- as.Date(pdt$Date, format=findDateFormat(pdt$Date))

ohc.breaks = seq(0, 1, length.out=201)
ohc.mid = ohc.breaks[1:(length(ohc.breaks)-1)]
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
ohc.col = jet.colors(length(ohc.breaks)-1) #[as.vector((dataT))]

# set axes
x.at <- pretty(xlims, 5)
x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
y.at <- pretty(ylims, 5)
y.labels <- parse(text=paste(abs(y.at), "*degree~N", sep=""))

load(paste(woa.dir,'woa.quarter.atl.rda',sep=''))

#for (i in 62:length(pdt.udates)){
  i = 50
  time <- pdt.udates[i]
  pdt.i <- pdt[which(pdt$Date == time),]
  
  # open day's hycom data
  nc <- RNetCDF::open.nc(paste(hycom.dir, ptt, '_', as.Date(time), '.nc', sep=''))
  dat <- RNetCDF::var.get.nc(nc, 'water_temp') * RNetCDF::att.get.nc(nc, 'water_temp', attribute='scale_factor') + 
    RNetCDF::att.get.nc(nc, variable='water_temp', attribute='add_offset')
  #if(i==1){
    lon <- RNetCDF::var.get.nc(nc, 'lon') - 360
    lat <- RNetCDF::var.get.nc(nc, 'lat')
    depth <- RNetCDF::var.get.nc(nc, 'depth')
  #}
  
  depIdx = unique(apply(as.data.frame(pdt.i$Depth), 1, FUN = function(x) which.min((x - depth) ^ 2)))
  hycomDep <- depth[depIdx]
  
  s <- raster::flip(raster::brick(dat[,,depIdx], xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), transpose=T), 2)
  
  # get locfit results
  # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
  suppressWarnings(fit.low <- locfit::locfit(pdt.i$MinTemp ~ pdt.i$Depth))
  suppressWarnings(fit.high <- locfit::locfit(pdt.i$MaxTemp ~ pdt.i$Depth))
  n = length(hycomDep)
  pred.low = stats::predict(fit.low, newdata = hycomDep, se = T, get.data = T)
  pred.high = stats::predict(fit.high, newdata = hycomDep, se = T, get.data = T)
  df = data.frame(low = pred.low$fit,# - pred.low$se.fit * sqrt(n),
                  high = pred.high$fit,# + pred.high$se.fit * sqrt(n),
                  depth = hycomDep)
  
  idx <- which(dateVec %in% pdt.udates[i])
  spot.prof <- raster::extract(s, cbind(spot$lon[idx], spot$lat[idx]))
  
  woa.lon.idx <- which.min(abs(spot$lon[idx] - woa.quarter$lon))
  woa.lat.idx <- which.min(abs(spot$lat[idx] - woa.quarter$lat))
  spot.woa <- woa.quarter$watertemp[woa.lon.idx, woa.lat.idx,,lubridate::month(time)]
  
  # start plotting
  #800 wide x 550 tall
  outDir <- '~/HMMoce_run/data/141254/check_hycom/'
  png(paste(outDir, ptt, '_prof_', dateVec[idx], '.png', sep=''),
      width=16, height=11, units='in', res=400)
  pdf(paste(outDir, ptt, '_prof_', dateVec[idx], '.pdf', sep=''),
      width=16, height=11)
  nf <- layout(matrix(c(1,2,7,
                        3,4,7,
                        5,6,7), 3, 3, byrow=T), widths=c(8,1.8,6), heights=c(5,5,5))
  #layout.show(nf)
  
  ohc.coords <- xyFromCell(L.ohc[[idx]], which.max(L.ohc[[idx]]))
  hyc.coords <- xyFromCell(L.hycom[[idx]], which.max(L.hycom[[idx]]))
  woa.coords <- xyFromCell(L.woa[[idx]], which.max(L.woa[[idx]]))
  woa.prof <- woa.quarter$watertemp[which.min(abs(woa.coords[1] - woa.quarter$lon)),
                                    which.min(abs(woa.coords[2] - woa.quarter$lat)),,
                                    lubridate::month(time)]
  
  # plot 1 - ohc map
  image(L.ohc[[idx]], maxpixels=ncell(L.ohc[[idx]]), axes=F, xlim=xlims, ylim=ylims,
        col=ohc.col, breaks=ohc.breaks, xlab='', ylab='')
  axis(1, at=x.at, labels=x.labels)
  axis(2, at=y.at, labels=y.labels)
  fields::world(add=T, fill=T, col='grey60')
  points(spot$lon[idx], spot$lat[idx], pch=23, col='black',bg='purple',cex=2.5)
  #points(ohc.coords, pch=21, col='white', bg=NULL, cex=2.5, lwd=2)
  text(ohc.coords, labels = 'A', col='white', cex=2, lwd=2)
  text(-83, 45, 'OHC')
  box()
  
  # plot 2 - ohc legend
  plot.new()
  #image(1, ohc.mid, t(as.matrix(ohc.mid)), breaks=ohc.breaks, col=ohc.col, axes=FALSE, xlab="",
  #      ylab='OHC likelihood')
  #axis(2);box();
  
  # plot 3 - hycom map
  image(L.hycom[[idx]], maxpixels=ncell(L.hycom[[idx]]), axes=F, xlim=xlims, ylim=ylims,
        col=ohc.col, breaks=ohc.breaks, xlab='', ylab='')
  axis(1, at=x.at, labels=x.labels)
  axis(2, at=y.at, labels=y.labels)
  fields::world(add=T, fill=T, col='grey60')
  #points(spot$lon[idx], spot$lat[idx], pch=23, col='black',bg='purple',cex=2.5)
  text(spot$lon[idx], spot$lat[idx], labels='B', col='white', cex=2, lwd=2)
  #points(hyc.coords, pch=21, col='white', bg=NULL, cex=2.5, lwd=2)
  text(hyc.coords[1], hyc.coords[2]+1, labels='C', col='white', cex=2, lwd=2)
  text(-83, 45, 'HYCOM')
  box()
  
  # plot 4 - hycom legend
  image(1, ohc.mid, t(as.matrix(ohc.mid)), breaks=ohc.breaks, col=ohc.col, axes=FALSE, xlab="",
        ylab='Likelihood')
  axis(2);box();
  
  # plot 5 - woa map
  image(L.woa[[idx]], maxpixels=ncell(L.woa[[idx]]), axes=F, xlim=xlims, ylim=ylims,
        col=ohc.col, breaks=ohc.breaks, xlab='', ylab='')
  axis(1, at=x.at, labels=x.labels)
  axis(2, at=y.at, labels=y.labels)
  fields::world(add=T, fill=T, col='grey60')
  #points(spot$lon[idx], spot$lat[idx], pch=23, col='black',bg='purple',cex=2.5)
  text(spot$lon[idx], spot$lat[idx], labels='D', col='white', cex=2, lwd=2)
  #points(woa.coords, pch=21, col='white', bg=NULL, cex=2.5, lwd=2)
  text(woa.coords, labels='E', col='white', cex=2, lwd=2)
  text(-83, 45, 'WOA')
  box()
  
  # plot 6 - woa legend
  plot.new()
  #image(1, ohc.mid, t(as.matrix(ohc.mid)), breaks=ohc.breaks, col=ohc.col, axes=FALSE, xlab="",
  #      ylab='WOA likelihood')
  #axis(2);box();
  
  # plot 7 - depth profiles
  maxZ <- max(pdt.i$Depth) + 50
  t.lims <- c(min(pdt.i$MinTemp) - 6.5, max(pdt.i$MaxTemp) + 1)
  t.at <- pretty(t.lims, 7)
  t.labels <- parse(text=paste(abs(t.at), "*degree~C", sep=""))
  
  ohc.prof <- raster::extract(s, cbind(xyFromCell(L.ohc[[idx]], which.max(L.ohc[[idx]]))))
  hyc.prof <- raster::extract(s, cbind(xyFromCell(L.hycom[[idx]], which.max(L.hycom[[idx]]))))
  if(nrow(ohc.prof) > 1) ohc.prof <- ohc.prof[1,]
  if(nrow(hyc.prof) > 1) hyc.prof <- hyc.prof[1,]
  
  yleg <- .9 * maxZ
  xleg <- .9 * t.lims[2]
  
  linecols <- RColorBrewer::brewer.pal(5,'Set1')
  
  plot(pdt.i$MinTemp, pdt.i$Depth, type='l', ylim=c(max(pdt.i$Depth), 0), xlim=t.lims,
       xlab='Water Temperature', ylab='Depth (m)', axes=F, col='white')
  axis(2)
  axis(1, at=t.at, labels=t.labels)
  box()
  px <- c(df$low, rev(df$high))
  py <- c(df$depth, rev(df$depth))
  polygon(px, py, col=adjustcolor('grey60', alpha.f=.7))
  #lines(pdt.i$MaxTemp, pdt.i$Depth)
  lines(spot.prof, hycomDep, col=linecols[1], lwd=1.5) # B
  lines(ohc.prof, hycomDep, col=linecols[2], lwd=1.5) # A
  lines(hyc.prof, hycomDep, col=linecols[3], lwd=1.5) # C
  lines(spot.woa[3:34], woa.quarter$depth[3:34], col=linecols[4], lwd=1.5) # D
  lines(woa.prof[3:34], woa.quarter$depth[3:34], col=linecols[5], lwd=1.5) # E
  text(ohc.prof[1], -2, 'A')
  text(spot.prof[1], -2, 'B')
  text(hyc.prof[1], -2, 'C')
  text(spot.woa[1], -2, 'D')
  text(woa.prof[1], -2, 'E')
  
  #lines(df$low, df$depth, lty=2, lwd=1.3)
  #lines(df$high, df$depth, lty=2, lwd=1.3)
  #legend(xleg, yleg,
  #       c('pdt','spot','ohcmax','hycmax','locfit'),
  #       lty=c(1,1,1,1,2),
  #       col=c('black','red','blue','green','black'))
  dev.off()
  
#}
  