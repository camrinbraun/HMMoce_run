# FIGURE 1 - BASKING SHARK MS
load('~/ebs/Data/BaskingSharks/batch/bask_results_v4.rda')
load('~/ebs/Data/BaskingSharks/batch/all_bask_tad_v3.rda')
library(raster); library(fields); library(rgdal)

# prepare seasonal RD data
# bask.res is a list of length 3
# [[1]] is named list (all ptts) with allRD, behavRD, and df of track for each
# [[2]] is raster stack of allRD for each ptt (each has own layer)
# [[3]] is df of all tracks together
seasonList <- list()
for (i in 1:length(bask.res[[1]])){
  seasonList[[i]] <- bask.res[[1]][[i]]$seasonRD
}

for (i in 1:length(bask.res[[1]])){
  # resample to large spatial limits
  seasonList[[i]] <- raster::resample(seasonList[[i]], seasonList[[29]])
}

allSeason <- list(stack(lapply(seasonList, FUN=function(x) x[[1]])),
                  stack(lapply(seasonList, FUN=function(x) x[[2]])),
                  stack(lapply(seasonList, FUN=function(x) x[[3]])),
                  stack(lapply(seasonList, FUN=function(x) x[[4]])))

for (i in 1:4){
  for (t in 1:37){
    allSeason[[i]][[t]] <- allSeason[[i]][[t]] / cellStats(allSeason[[i]][[t]], 'max')
  }
}

allSeason <- stack(lapply(allSeason, FUN=function(x) sum(x, na.rm=T)))

# prepare TAD summarized data
summ <- all.tad.res$summ

# set colors
require(RColorBrewer)
chl.pal <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))#(100))
zmin <- min(cellStats(allSeason, 'min'), na.rm=T)
zmax <- max(cellStats(allSeason, 'max'), na.rm=T)
zbreaks <- seq(0.01, 1, length.out=201)#, zmax)
zmid <- zbreaks[1:(length(zbreaks)-1)]
zcol <- chl.pal(length(zbreaks)-1)

# set lims
xlims=c(-85,-20); ylims=c(-20,50)
x.at <- pretty(xlims, 5)
x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
y.at <- pretty(ylims, 5)
y.labels <- c(parse(text=paste(abs(y.at[which(y.at < 0)]), "*degree~S", sep="")), parse(text=paste(y.at[which(y.at >= 0)], "*degree~N", sep="")))

old.mar <- par()$mar
mar.default <- par('mar')
mar.r1 <- mar.default
mar.r2 <- mar.default
mar.r1 <- c(0,4,4,0)
mar.r2 <- c(4,4,4,0)

# BUILD THE PLOT
pdf('~/ebs/Data/BaskingSharks/batch/fig2_bask_season.pdf', width=16, height=10)
nf <- layout(matrix(c(1,2,3,4,
                      5,6,7,8), nrow=2, ncol=4, byrow=T), widths=c(5,3,5,3), heights=c(5,5))
#layout.show(nf)
#par(fig=c(0,1,0,1))
#par(mar=mar.left, fig=c(0,.59,0,1), new=T)

# spring - plot 1
par(mar=mar.r1)
image(allSeason[[1]] / cellStats(allSeason[[1]], 'max'), maxpixels=ncell(allSeason[[1]]), xlim=xlims, ylim=ylims,
      col=zcol, breaks=zbreaks, xlab='', ylab='', axes=F)
contour(allSeason[[1]] / cellStats(allSeason[[1]], 'max'), add=T, levels=c(.25, .5), labels=c('75%','50%'), lty=c(2,1))
fields::world(add=T, fill=T, col='grey80', border='grey80')
axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels)
text(-22,-18,'A', font=2)
box()

# spring - plot 2
par(mar=c(0,4,4,2))
b <- barplot(rev(summ[which(summ$season==1),3]), horiz=T, axes=F, xlim=c(0,67))
axis(2, at=b, labels=rev(c('<10','10-25','25-50','50-200','200-400','400-1000','>2000')))
axis(1)
text(59,0.17,'B', font=2)
#mtext('Time-at-depth (%)', side=1, line=2.7)
mtext('Depth (m)', side=2, line=2.25, cex=.8)

# summer - plot 3
par(mar=mar.r1)
image(allSeason[[2]] / cellStats(allSeason[[2]], 'max'), maxpixels=ncell(allSeason[[1]]), xlim=xlims, ylim=ylims,
      col=zcol, breaks=zbreaks, xlab='', ylab='', axes=F)
contour(allSeason[[2]] / cellStats(allSeason[[2]], 'max'), add=T, levels=c(.25, .5), labels=c('75%','50%'), lty=c(2,1))
fields::world(add=T, fill=T, col='grey80', border='grey80')
axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels)
text(-22,-18,'C', font=2)
box()

# summer - plot 4
par(mar=c(0,4,4,2))
b <- barplot(rev(summ[which(summ$season==2),3]), horiz=T, axes=F, xlim=c(0,67))
axis(2, at=b, labels=rev(c('<10','10-25','25-50','50-200','200-400','400-1000','>2000')))
axis(1)
text(59,0.17,'D', font=2)
#mtext('Time-at-depth (%)', side=1, line=2.7)
mtext('Depth (m)', side=2, line=2.25, cex=.8)

# fall - plot 5
par(mar=mar.r2)
image(allSeason[[3]] / cellStats(allSeason[[3]], 'max'), maxpixels=ncell(allSeason[[1]]), xlim=xlims, ylim=ylims,
      col=zcol, breaks=zbreaks, xlab='', ylab='', axes=F)
contour(allSeason[[3]] / cellStats(allSeason[[3]], 'max'), add=T, levels=c(.25, .5), labels=c('75%','50%'), lty=c(2,1))
fields::world(add=T, fill=T, col='grey80', border='grey80')
axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels)
text(-22,-18,'E', font=2)
box()

# fall - plot 6
par(mar=c(4,4,2,2))
b <- barplot(rev(summ[which(summ$season==3),3]), horiz=T, axes=F, xlim=c(0,67))
axis(2, at=b, labels=rev(c('<10','10-25','25-50','50-200','200-400','400-1000','>2000')))
axis(1)
mtext('Time-at-depth (%)', side=1, line=2.5, cex=.8)
mtext('Depth (m)', side=2, line=2.25, cex=.8)
text(59,0.17,'F', font=2)

# winter - plot 7
par(mar=mar.r2)
image(allSeason[[4]] / cellStats(allSeason[[4]], 'max'), maxpixels=ncell(allSeason[[1]]), xlim=xlims, ylim=ylims,
      col=zcol, breaks=zbreaks, xlab='', ylab='', axes=F)
contour(allSeason[[4]] / cellStats(allSeason[[4]], 'max'), add=T, levels=c(.25, .5), labels=c('75%','50%'), lty=c(2,1))
fields::world(add=T, fill=T, col='grey80', border='grey80')
axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels)
text(-22,-18,'G', font=2)
box()
  
# winter - plot 8
par(mar=c(4,4,2,2))
b <- barplot(rev(summ[which(summ$season==4),3]), horiz=T, axes=F, xlim=c(0,67))
axis(2, at=b, labels=rev(c('<10','10-25','25-50','50-200','200-400','400-1000','>2000')))
axis(1)
mtext('Time-at-depth (%)', side=1, line=2.5, cex=.8)
mtext('Depth (m)', side=2, line=2.25, cex=.8)
text(59,0.17,'H', font=2)

dev.off()
  
  
  