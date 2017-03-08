
require(raster); require(ggplot2); require(gridExtra); require(rgdal);
require(maptools); require(plyr); require(HMMoce); require(fields);
require(maps)

source('~/Documents/WHOI/RCode/hmm_ms/R/plot_funs.r')
load('~/Documents/WHOI/RData/Blues/2015/141259/blue259_allL.RData')

# PROVIDE FIXED KERNEL PARAMETERS
par0 <- c(10,10,3,1)
D1 <- par0[1:2] # parameters for kernel 1. this is migratory behavior mode
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode

# GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
K1 <- gausskern(D1[1], D1[2], muadv = 0)
K2 <- gausskern(D2[1], D2[2], muadv = 0)

#-------------------------#
# prep the rasters
#-------------------------#
i = 101 #"2016-01-27"
r1 <- L.res[[1]]$L.light[[i]]
r2 <- L.res[[1]]$L.sst[[i]]
r3 <- L.res[[1]]$L.ohc[[i]]
r4 <- L.res[[1]]$L.prof[[12]]
plot(r4)
r4 <- r4 + 2
r4 <- r4 / cellStats(r4, 'max')
r4[r4 < .96] <- NA
r4;plot(r4)
# convert the rasters to points for plotting
map.p <- rasterToPoints(r1)
#Make the points a dataframe for ggplot
df <- data.frame(map.p)
#Make appropriate column headings
colnames(df) <- c('Longitude', 'Latitude', 'MAP')
df$group = 1
df1 <- df
df1$MAP <- df1$MAP / max(df1$MAP, na.rm=T)
df1$MAP[df1$MAP < .01] <- NA

# convert the rasters to points for plotting
map.p <- rasterToPoints(r2)
#Make the points a dataframe for ggplot
df <- data.frame(map.p)
#Make appropriate column headings
colnames(df) <- c('Longitude', 'Latitude', 'MAP')
df$group = 1
df2 <- df
df2$MAP <- df2$MAP / max(df2$MAP, na.rm=T)
df2$MAP[df2$MAP < .01] <- NA

# convert the rasters to points for plotting
map.p <- rasterToPoints(r3)
#Make the points a dataframe for ggplot
df <- data.frame(map.p)
#Make appropriate column headings
colnames(df) <- c('Longitude', 'Latitude', 'MAP')
df$group = 1
df3 <- df
df3$MAP <- df3$MAP / max(df3$MAP, na.rm=T)
df3$MAP[df3$MAP < .01] <- NA

# convert the rasters to points for plotting
map.p <- rasterToPoints(r4)
#Make the points a dataframe for ggplot
df <- data.frame(map.p)
#Make appropriate column headings
colnames(df) <- c('Longitude', 'Latitude', 'MAP')
df$group = 1
df4 <- df
df4$MAP <- df4$MAP / max(df4$MAP, na.rm=T)
df4$MAP[df4$MAP < .01] <- NA

#-------------------------#
# prep the countries
#-------------------------#
#countries <- readShapePoly('/Users/Cam/Documents/WHOI/RData/Countries/countries_dissolve.shp')
#proj4string(countries)=CRS("+init=epsg:3395") #world mercator projection
setwd('~/Documents/WHOI/RData/Countries/')
countries = readOGR(dsn=".", layer="countries", verbose = FALSE)
countries@data$id = rownames(countries@data)
countries.points = fortify(countries, region="id")
countries.df = join(countries.points, countries@data, by="id")
#setwd(wd)

#-------------------------#
# make the plots
#-------------------------#
lims <- c(-81,-45,24,44)

p1 <- ggplot(countries.df) + aes(long,lat,group=group) + geom_polygon() +
  geom_raster(data=df1, aes(y=Latitude, x=Longitude,fill=MAP)) +
  scale_fill_gradientn(colours = jet.colors(30), na.value='white') + xlab('') + ylab('') +
  #geom_point(-73.4393, 28.4547) +
  coord_equal() + coord_cartesian(xlim=c(lims[1],lims[2]), ylim=c(lims[3],lims[4])) +
  theme_bw() + theme(panel.border = element_rect(colour='black',size=1.15), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     legend.position='none',axis.text.x=element_blank())

p2 <- ggplot(countries.df) + aes(long,lat,group=group) + geom_polygon() +
  geom_raster(data=df2, aes(y=Latitude, x=Longitude,fill=MAP)) +
  scale_fill_gradientn(colours = jet.colors(30), na.value='white') + xlab('') + ylab('') +
  coord_equal() + coord_cartesian(xlim=c(lims[1],lims[2]), ylim=c(lims[3],lims[4])) +
  theme_bw() + theme(panel.border = element_rect(colour='black',size=1.15), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     axis.text.x=element_blank(),legend.title=element_blank(),
                     axis.text.y=element_blank())

p3 <- ggplot(countries.df) + aes(long,lat,group=group) + geom_polygon() +
  geom_raster(data=df3, aes(y=Latitude, x=Longitude,fill=MAP)) +
  scale_fill_gradientn(colours = jet.colors(30), na.value='white') + xlab('') + ylab('') +
  coord_equal() + coord_cartesian(xlim=c(lims[1],lims[2]), ylim=c(lims[3],lims[4])) +
  theme_bw() + theme(panel.border = element_rect(colour='black',size=1.15), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     legend.position='none')

p4 <- ggplot(countries.df) + aes(long,lat,group=group) + geom_polygon() +
  geom_raster(data=df4, aes(y=Latitude, x=Longitude,fill=MAP)) +
  scale_fill_gradientn(colours = jet.colors(30), na.value='white') + xlab('') + ylab('') +
  coord_equal() + coord_cartesian(xlim=c(lims[1],lims[2]), ylim=c(lims[3],lims[4])) +
  theme_bw() + theme(panel.border = element_rect(colour='black',size=1.15), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_blank(),
                     legend.position='none', axis.text.y=element_blank())

legend <- g_legend(p2)

grid.arrange(p1, p2+theme(legend.position='none'), p3, p4, legend, ncol = 3,
             layout_matrix = cbind(c(1,3), c(2,4), c(5)),
             widths = c(3,3,.85))
