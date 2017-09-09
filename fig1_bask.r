# FIGURE 1 - BASKING SHARK MS
load('~/ebs/Data/BaskingSharks/batch/bask_results_v4.rda')
library(raster); library(fields); library(rgdal)

# bask.res is a list of length 3
# [[1]] is named list (all ptts) with allRD, behavRD, and df of track for each
# [[2]] is raster stack of allRD for each ptt (each has own layer)
# [[3]] is df of all tracks together

# dlist is density list of lat by month
dlist <- list()
for(i in c(1:12)){
  dlist[[i]] <- density(locs$lat[which(locs$month == i)], 
                          bw=1.5, from=min(locs$lat, na.rm=T),
                          to=max(locs$lat, na.rm=T))
  dlist[[i]]$x <- c(ylims[1], dlist[[i]]$x, ylims[2]+1)
  dlist[[i]]$y <- c(0, dlist[[i]]$y, 0)

  dlist[[i]]$n.indiv <- length(unique(locs$ptt[which(locs$month == i)]))
  dlist[[i]]$month <- i
}


locs <- bask.res[[3]]
month.colors=c('#3B60AD','#43B883','#D8E47E','#FCCF52','#F47A2D','#E61E25','#F04F2E','#FCB345','#F5EE6C','#C0D48A','#82B3A9','#2A93C7')
locs$month <- lubridate::month(locs$date)
locs$fill <- month.colors[locs$month]

xlims=c(-85,-25); ylims=c(-10,45)
x.at <- pretty(xlims, 5)
x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
y.at <- pretty(ylims, 5)
y.labels <- c(parse(text=paste(abs(y.at[which(y.at < 0)]), "*degree~S", sep="")), parse(text=paste(y.at[which(y.at >= 0)], "*degree~N", sep="")))

old.mar <- par()$mar
mar.default <- par('mar')
mar.left <- mar.default
mar.right <- mar.default
mar.left[4] <- 0
mar.right[2] <- 0

# BUILD THE PLOT
pdf('~/ebs/Data/BaskingSharks/batch/fig1_bask_v4.pdf', width=12, height=12)
#nf <- layout(matrix(c(1,2,
#                      1,2), 1, 2, byrow=T), widths=c(5,3), heights=c(5))
#layout.show(nf)
par(fig=c(0,1,0,1))
par(mar=mar.left, fig=c(0,.59,0,1), new=T)

plot(locs$lon, locs$lat, type='n', xlim=xlims, ylim=c(ylims[1],ylims[2]+1), axes=F, xlab='', ylab='')
#lines(gs@lines[[1]]@Lines[[10]],lty=3) #gulf stream contour
#lines(c200, lty=2) # 200m contour
fields::world(add=T, fill=T, col='grey80', border='grey80')

plot.lines <- split(locs, locs$ptt) # add track lines
lapply(plot.lines, function(z) lines(z$lon, z$lat))

#na.idx <- which(is.na(locs$fill))
#points(locs$lon[na.idx],locs$lat[na.idx], pch=21, bg='grey80', col='grey80', cex=1)
points(locs$lon, locs$lat, pch=21, bg=locs$fill, col=locs$fill, cex=1)

lapply(plot.lines, function(z) points(z$lon[1], z$lat[1], pch=24, bg='green', cex=1.5))
lapply(plot.lines, function(z) points(z$lon[nrow(z)], z$lat[nrow(z)], pch=25, bg='red', cex=1.5))

img <- readPNG("~/ebs/Data/BaskingSharks/batch/month_point_legend_r-01.png")
rasterImage(img,-85,-11,-75,9)

axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels);
box()

# panel 2, density of lat by month
par(mar=mar.right, fig=c(.60,1,0,1), new=T)
plot(dlist[[1]]$y, dlist[[1]]$x, type='n', axes=F, xlab='', ylab='', xlim=c(0,.35), ylim=c(ylims[1], ylims[2]+1))

for(i in c(12:1)){
  polygon(dlist[[i]]$y+((i-1)*.022), dlist[[i]]$x, col=rgb(204/255,204/255,204/255,.7), border=NULL)
  #lines(dlist[[i]]$y, dlist[[i]]$x)
}
xvec <- seq(0, .022*11, by=.022)
#xvec <- xvec[-8]
yvec <- rep(47, length.out=12)
yvec2 <- rep(-11.25, length.out=12)
labels <- c('J','F','M','A','M','J','J','A','S','O','N','D')
n.month <- unlist(lapply(dlist, FUN=function(x) x$n.indiv))
text(xvec, yvec, labels)
text(xvec, yvec2, n.month)
#mtext('# of individuals', 1, at=-45)
#mtext('Month', 3, at=.75)
#mtext('Month', side=3, adj=.5)

dev.off()


#=====================
