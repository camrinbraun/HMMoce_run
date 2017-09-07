# fig 3 is two panel plot of bask tracks
# panel A has ptts: 88137, 88141, 88146
listA <- c(88137, 88141, 88146)
# panel b has: 52562, 100975, 110502
listB <- c(52562, 100975, 110502)
# bathy underneath both with bath.colors, grey landmass, tag/pop as triangles?

library(raster); library(fields); library(RColorBrewer)
#bathy <- raster('~/ebs/EnvData/bathy/BaskingSharks/bask_bathy_big.gri')
bathy <- HMMoce::get.bath.data(-100, -10, -30, 60, res = c(.5), seaonly=F)
writeRaster(bathy, '~/ebs/EnvData/bathy/BaskingSharks/plot_big_bathy.grd')
bathy <- raster('~/ebs/EnvData/bathy/BaskingSharks/plot_big_bathy.grd')

load('~/ebs/Data/BaskingSharks/batch/bask_results_v2.rda')
df <- bask.res[[3]]; rm(bask.res)
df$month <- lubridate::month(df$date)
df.a <- df[which(df$ptt %in% listA),]
df.b <- df[which(df$ptt %in% listB),]

# set colors
bath.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan"))
month.colors=c('#3B60AD','#43B883','#D8E47E','#FCCF52','#F47A2D','#E61E25','#F04F2E','#FCB345','#F5EE6C','#C0D48A','#82B3A9','#2A93C7')

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

# plot 1
# set lims
xlims=c(-85,-50); ylims=c(10,47)
x.at <- pretty(xlims, 5)
x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
y.at <- pretty(ylims, 5)
y.labels <- c(parse(text=paste(abs(y.at[which(y.at < 0)]), "*degree~S", sep="")), parse(text=paste(y.at[which(y.at >= 0)], "*degree~N", sep="")))

par(mar=mar.r1)
image(bathy, maxpixels=ncell(bathy), xlim=xlims, ylim=ylims,
      col=bath.colors(255), xlab='', ylab='', axes=F)
fields::world(add=T, fill=T, col='grey80', border='grey80')
axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels)
text(-22,-18,'A', font=2)
box()

points(df.a$lon, df.a$lat, pch=21, bg=month.colors[df.a$month])

# set lims
xlims=c(-85,-20); ylims=c(-20,50)
x.at <- pretty(xlims, 5)
x.labels <- parse(text=paste(abs(x.at), "*degree~W", sep=""))
y.at <- pretty(ylims, 5)
y.labels <- c(parse(text=paste(abs(y.at[which(y.at < 0)]), "*degree~S", sep="")), parse(text=paste(y.at[which(y.at >= 0)], "*degree~N", sep="")))

image(bathy, maxpixels=ncell(bathy), xlim=xlims, ylim=ylims,
      col=bath.colors(255), xlab='', ylab='', axes=F)
fields::world(add=T, fill=T, col='grey80', border='grey80')
axis(1, at=x.at, labels=x.labels)
axis(2, at=y.at, labels=y.labels)
text(-22,-18,'A', font=2)
box()

points(df.b$lon, df.b$lat, pch=21, bg=month.colors[df.b$month])



