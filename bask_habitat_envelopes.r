#=================================================================#
# BASK HABITAT ENVELOPES
#=================================================================#
library(plyr); library(ggplot2)

sar <- extent(c(xmn=-78, xmx=-50, ymn=25, ymx=38))
ne <- extent(c(xmn=-73, xmx=-65, ymn=40, ymx=44))
ant <- extent(c(xmn=-70, xmx=-55, ymn=10, ymx=25))
sam <- extent(c(xmn=-60, xmx=-20, ymn=-12, ymx=10))

plot(df$lon, df$lat, xlim=c(-75,-60), ylim=c(35,50))
world(add=T)
plot(sar, add=T, col='green')
plot(ne, add=T, col='yellow')
plot(ant, add=T, col='blue')
plot(sam, add=T, col='red')

df$area <- NA
df$area[which(df$lon <= sar[2] & df$lon >= sar[1] &
                df$lat >= sar[3] & df$lat <= sar[4])] <- 'SARG'
df$area[which(df$lon <= ne[2] & df$lon >= ne[1] &
                df$lat >= ne[3] & df$lat <= ne[4])] <- 'NENG'
df$area[which(df$lon <= ant[2] & df$lon >= ant[1] &
                df$lat >= ant[3] & df$lat <= ant[4])] <- 'ANTL'
df$area[which(df$lon <= sam[2] & df$lon >= sam[1] &
                df$lat >= sam[3] & df$lat <= sam[4])] <- 'SAMR'
df$area <- as.factor(df$area)

for (ii in 2004:2012){
  # WINTER
  if(ii != 2004) df$season[which(df$date >= paste(ii-1,'-12-20',sep='') & df$date < paste(ii-1,'-12-31',sep=''))] <- 4
  df$season[which(df$date >= paste(ii,'-01-01',sep='') & df$date < paste(ii,'-03-20',sep=''))] <- 4
  
  # FALL
  df$season[which(df$date >= paste(ii,'-09-20',sep='') & df$date < paste(ii,'-12-20',sep=''))] <- 3
  
  #SPR
  df$season[which(df$date >= paste(ii,'-03-20',sep='') & df$date < paste(ii,'-06-20',sep=''))] <- 1
  
  #SUM
  df$season[which(df$date >= paste(ii,'-06-20',sep='') & df$date < paste(ii,'-09-20',sep=''))] <- 2
  
}
df$season[which(is.na(df$season))] <- 4

allpdt$season[which(allpdt$season == 1)] <- 'Spring'
allpdt$season[which(allpdt$season == 2)] <- 'Summer'
allpdt$season[which(allpdt$season == 3)] <- 'Fall'
allpdt$season[which(allpdt$season == 4)] <- 'Winter'


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
  pdt <- read.wc(ptt, filename = paste(myDir, ptt, '-PDTs.csv',sep=''), type = 'pdt', tag=tag, pop=pop); 
  pdt.udates <- pdt$udates; pdt <- pdt$data
  tail(pdt)
  print(meta$PopDate[ii])
  pdt$Date <- as.Date(pdt$Date, format = findDateFormat(pdt$Date))
  names(pdt) <- tolower(names(pdt))
  keepNames <- c('ptt','date','binnum','depth','mintemp','maxtemp')
  pdt <- pdt[,which(names(pdt) %in% keepNames)]
  
  df.i <- df[which(df$ptt == ptt),]
  
  pdt.i <- merge(pdt, df.i, by='date')
  
  if(ii==1){
    allpdt <- pdt.i
  } else{
    allpdt <- rbind(allpdt, pdt.i)
  }
  
}

allpdt.save <- allpdt
allpdt <- allpdt[which(!is.na(allpdt$area)),]
allpdt$depth <- allpdt$depth * -1

allpdt.dat <- within(allpdt, area <- factor(area, levels = c('NENG','SARG','ANTL','SAMR')))
allpdt.dat <- within(allpdt.dat, season <- factor(season, levels = c('Fall','Winter','Spring','Summer')))


jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF",
                                 "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

plot_labeller <- function(variable,value){
  if (variable == 'season') {
    value[value=="1"] <- "Spring"
    value[value=="2"] <- "Summer"
    value[value=="3"] <- "Fall"
    value[value=="4"] <- "Winter"
    } else {
    value[value=='ANTL'] <- 'Antilles'
    value[value=='NENG'] <- 'New England'
    value[value=='SAMR'] <- 'S. America'
    value[value=='SARG'] <- 'Sargasso'
  }
}

mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var == "season") { 
    value[value=="1"] <- "Spring"
    value[value=="2"] <- "Summer"
    value[value=="3"] <- "Fall"
    value[value=="4"] <- "Winter"
  }
  return(value)
}

allpdt.dat$midtemp <- (allpdt.dat$maxtemp + allpdt.dat$mintemp) / 2

ylab <- 'Depth (m)'
xlab <- expression(paste("Temperature (",degree,"C)"))
#  bzz <- c(1,5,10,20,50,100,200,500,1000,2000,4000,6000)
bzz <- c(1,5,10,20,50,100,200)
byy <- rev(seq(0,1000,100)*-1)
#ly <- rep("",length(byy))
ly <- as.character(byy*-1)
#ly[seq(1,length(byy),by=100)] <- byy[seq(1,length(byy),by=100)]*-1
ly[1] <- ">1000"
bxx <- seq(4,32,2)
lx <- rep("",length(bxx))
lx[seq(2,length(bxx),by=2)] <- bxx[seq(2,length(bxx),by=2)]
sx <- seq(4, 32, by = 1)
sy <- rev(seq(0, 1000, by = 25)*-1)
p <- ggplot(allpdt.dat, aes(x = midtemp, y = depth))
p <- p + geom_bin2d(breaks = list(x = sx, y = sy)) 
p <- p + scale_x_continuous(name=xlab, breaks= bxx, labels=format(lx)) 
p <- p + scale_y_continuous(name=ylab, breaks= byy, labels=format(ly))
p <- p + scale_fill_gradientn(" Frequency of\n occurence\n", trans = 'log', colours = jet.colors(10), 
                              breaks = bzz, labels = format(bzz), 
                              guide = guide_colorbar(barwidth = 1.5, barheight = 30))
#p + facet_grid(. ~ cntry, labeller=mf_labeller) + theme_bw(base_size = 10)
p <- p + facet_grid(season ~ area) + theme_bw(base_size = 10) + theme(panel.grid=element_blank())
atxt <- ddply(.data=allpdt.dat, 
              .(season, area), 
              summarize, 
              n=paste("n =", length(midtemp)))
p + geom_text(data=atxt, aes(x=27.5, y=-1000, label=n),colour="black", size = 2.7, inherit.aes=FALSE, parse=FALSE)
dev.copy2pdf(width=10,height=10, file="~/ebs/Data/BaskingSharks/batch/HabitatEnvelopeSeasonal.pdf")

save.image('~/ebs/Data/BaskingSharks/batch/bask_hab_envelopes.rda')

# some descriptive stats
which(allpdt$area == 'NENG')


## END
#===========================
dat <- hab
ss <- read.csv("C:\\Docs\\StripedMarlin\\seasons.csv")
ss$stamp <- with(ss, paste(loc,month,sep="-"))
dat$season <- with(dat, paste(loc,month,sep="-"))
dat$season <- ss$season[match(dat$season,ss$stamp)]

{
  names(dat) <- tolower(names(dat))
  dat$depth[which(dat$depth<0)] = 0
  # Some NAs
  dat <- dat[-unique(c(which(is.na(dat$depth)),which(is.na(dat$temp)))),]
  dat$Depth <- dat$depth*-1
  dat$Depth[which(dat$Depth>360)] = 360
  
  country.code <- function(x,pos=1){
    x$cntry <- NA
    x$cntry[which(x[,pos]=="Australia")] <- "1"
    x$cntry[which(x[,pos]=="California")] <- "4"
    x$cntry[which(x[,pos]=="Hawaii")] <- "3"
    x$cntry[which(x[,pos]=="Mexico")] <- "5"
    x$cntry[which(x[,pos]=="New Zealand")] <- "2"
    x$cntry[which(x[,pos]=="Costa Rica")] <- "7"
    x$cntry[which(x[,pos]=="Ecuador")] <- "6"
    x$cntry[which(x[,pos]=="Panama")] <- "8" 
    return(x)
  }
  dat <- country.code(dat)
  
  dat$season <- factor(dat$season, levels = c("Summer", "Autumn", "Winter", "Spring"))
  
  mf_labeller <- function(var, value){
    value <- as.character(value)
    if (var=="cntry") { 
      value[value=="1"] <- "Australia"
      value[value=="4"] <- "California"
      value[value=="3"] <- "Hawaii"
      value[value=="5"] <- "Mexico"
      value[value=="2"] <- "New Zealand"
      value[value=="7"] <- "Costa Rica"
      value[value=="6"] <- "Ecuador"
      value[value=="8"] <- "Panama"		
    }
    return(value)
  }
  
  p <- ggplot(dat, aes_string(x = "temp", y = "Depth"))
  ylab <- 'Depth (m)'
  xlab <- expression(paste("Temperature (",degree,"C)"))
  #  bzz <- c(1,5,10,20,50,100,200,500,1000,2000,4000,6000)
  bzz <- c(1,5,10,20,50,100,200,500,1000,2000,3000)
  byy <- rev(seq(0,350,10)*-1)
  ly <- rep("",length(byy))
  ly[seq(1,length(byy),by=5)] <- byy[seq(1,length(byy),by=5)]*-1
  ly[1] <- ">350"
  bxx <- seq(8,32,2)
  lx <- rep("",length(bxx))
  lx[seq(2,length(bxx),by=2)] <- bxx[seq(2,length(bxx),by=2)]
  sx <- seq(9, 32, by = 1)
  sy <- rev(seq(0, 360, by = 10)*-1)
  p <- p + stat_bin2d(breaks = list(x = sx, y = sy)) 
  p <- p + scale_x_continuous(name=xlab, breaks= bxx, labels=format(lx)) 
  p <- p + scale_y_continuous(name=ylab, breaks= byy, labels=format(ly))
  p <- p + scale_fill_gradientn(" Frequency of\n occurence\n", trans = 'log', colours = jet.colors(10), 
                                breaks = bzz, labels = format(bzz), 
                                guide = guide_colorbar(barwidth = 1.5, barheight = 30))
  #p + facet_grid(. ~ cntry, labeller=mf_labeller) + theme_bw(base_size = 10)
  p <- p + facet_grid(season ~ cntry, labeller=mf_labeller) + theme_bw(base_size = 10)
  atxt <- ddply(.data=dat, 
                .(season, cntry), 
                summarize, 
                n=paste("n =", length(temp)))
  p + geom_text(data=atxt, aes(x=27.5, y=-360, label=n),colour="black", size = 2.7, inherit.aes=FALSE, parse=FALSE)
  dev.copy2pdf(width=12,height=10, file="C:\\Docs\\StripedMarlin\\figures-potential\\MlsHabitatEnvelopeSeasonal.pdf")
}
