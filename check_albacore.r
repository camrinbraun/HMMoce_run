## 1090251

work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
setwd(work_dir)


## check hycom at tagging location
hy <- raster::brick(paste0('./tmp/hycom/hycom_', tag, '.nc'))
hy <- raster::flip(hy, 'y')
prof_start <- raster::extract(hy, cbind(iniloc$lon[1], iniloc$lat[1]))
prof_end <- raster::extract(hy, cbind(iniloc$lon[2], iniloc$lat[2]))

plot(flip(hy[[1]],'y'))
world(add=T, wrap=c(0,360))
points(cbind(iniloc$lon[1], iniloc$lat[1]))
points(cbind(iniloc$lon[2], iniloc$lat[2]))

data_dir <- paste0('~/ebs/Data/albacore/', meta$instrument_name[i], '/cdb/')
etuff_file <- paste(data_dir, meta$instrument_name[i], '_eTUFF.txt', sep='')

etuff <- read_archival(etuff_file)
series <- get_series(etuff)

## depth-temp profiles: Date, Depth, MeanTemp
pdt <- data.frame(series %>% filter(!is.na(temperature)))
pdt$Date <- as.Date(pdt$DateTime)
hycom_depth <- c(0, 2, 4, 6, 8, 10, 12, 15, 20, 25,
                 30, 35, 40, 45, 50, 60, 70, 80, 90,
                 100, 125, 150, 200, 250, 300, 350, 
                 400, 500, 600, 700, 800, 900, 1000,
                 1250, 1500, 2000, 2500, 3000, 4000, 5000)
pdt$depth[which(pdt$depth < 0)] <- 0
pdt$depthInterval <- findInterval(pdt$depth, hycom_depth)
pdt <- data.frame(pdt %>% group_by(Date, depthInterval) %>% 
                    summarise(MinTemp = min(temperature), MaxTemp = max(temperature), 
                              MeanTemp = mean(temperature), n=n()))
pdt$depth <- hycom_depth[pdt$depthInterval]
pdt <- pdt[,c('Date','depth','MinTemp','MaxTemp')]
names(pdt)[2] <- 'Depth'
pdt$MinTemp <- round(pdt$MinTemp, 1)
pdt$MaxTemp <- round(pdt$MaxTemp, 1)
pdt$Date <- as.POSIXct(pdt$Date, tz='UTC')

idx <- which(pdt$Date == dateVec[2])
plot(c(prof_start), hycom_depth, ylim=c(100,0))
lines(pdt$MinTemp[idx], pdt$Depth[idx], col='blue')
lines(pdt$MaxTemp[idx], pdt$Depth[idx], col='red')

pdt.i <- pdt[idx,]
pdt.i$useTemp <- (pdt.i$MaxTemp + pdt.i$MinTemp) / 2

names(pdt.i) <- tolower(names(pdt.i))
suppressWarnings(
  fit.low <- locfit::locfit(pdt.i$mintemp ~ pdt.i$depth, maxk=500)
)
suppressWarnings(
  fit.high <- locfit::locfit(pdt.i$maxtemp ~ pdt.i$depth, maxk=500)
)
n = length(hycomDep)

#suppressWarnings(
pred.low = stats::predict(fit.low, newdata = hycomDep, se = TRUE, get.data = TRUE)
#suppressWarnings(
pred.high = stats::predict(fit.high, newdata = hycomDep, se = TRUE, get.data = TRUE)

df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                high = pred.high$fit + pred.high$se.fit * sqrt(n),
                depth = hycomDep)

minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = TRUE) / 10000
maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = TRUE) / 10000

dat <- as.array(hy)
dat <- dat[,,depIdx] - isotherm
ohc <- cp * rho * apply(dat, 1:2, sum, na.rm = TRUE) / 10000 
ohc[ohc == 0] <- NA

nc1 <- RNetCDF::open.nc(paste0('./tmp/hycom/hycom_', tag, '.nc'))
lon <- RNetCDF::var.get.nc(nc1, 'longitude')
lat <- RNetCDF::var.get.nc(nc1, 'latitude')

crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
if(!any(lon < 180)) lon <- lon - 360
list.ohc <- list(x = lon, y = lat, z = ohc)
ex <- raster::extent(list.ohc)
ohc <- raster::raster(list.ohc$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], crs)
#ohc <- raster::flip(ohc, direction = 'y')
plot(ohc)


## CONCLUSION -> all calcs appear to be working correctly. HYCOM is just off at depth in that area causing a shift in OHC way to the south

#===============
## try L.sst with WC-derived SST rather than eTUFF

sst <- data.frame(series %>% filter(!is.na(temperature)) %>% 
                    group_by(as.Date(DateTime)) %>%
                    summarise(n=n(),
                              length_temp = length(temperature[which.min(depth)]),
                              which_min = which.min(depth),
                              sst = temperature[which.min(depth)],
                              sst_depth = depth[which.min(depth)],
                              sst_maxT = max(temperature, na.rm=T),
                              sst_mean = mean(temperature[which(depth <= 2)])))
sst <- sst %>% dplyr::select(as.Date.DateTime., sst, sst_depth)
names(sst) <- c('Date','Temperature', 'sst_depth')
sst$Date <- as.POSIXct(sst$Date, tz='UTC')

wc_sst_dir <- paste0('~/ebs/Data/albacore/', meta$instrument_name[i], '/cdb/')
fList <- list.files(wc_sst_dir, full.names = TRUE)
wc_sst <- data.table::fread(fList[grep('-SST.csv', fList)], sep=',')
wc_sst$Date <- as.POSIXct(wc_sst$Date, format='%H:%M:%S %d-%b-%Y', tz='UTC')
wc_sst <- wc_sst %>% group_by(as.Date(Date)) %>% summarise(n=n(), mean_sst = mean(Temperature), min_sst = min(Temperature), max_sst = max(Temperature))
wc_sst$Date <- as.POSIXct(wc_sst$`as.Date(Date)`, tz='UTC')

#sst <- sst %>% group_by(as.Date(Date)) %>% summarise(n=n(), mean_sst = mean(Temperature), min_sst = min(Temperature), max_sst = max(Temperature))
sst2 <- archival_to_etuff(etuff$etuff, vars ='sst')
sst2 <- sst2 %>% group_by(as.Date(DateTime)) %>% summarise(n=n(), Temperature = mean(sst, na.rm=T), mean_sst = mean(sstMean, na.rm=T))
sst2$Date <- as.POSIXct(sst2$`as.Date(DateTime)`, tz='UTC')
sst2 <- data.frame(sst2)
for (b in 1:nrow(sst2)){
  if(is.nan(sst2$Temperature[b]) & !is.nan(sst2$mean_sst[b])) sst2$Temperature[b] <- sst2$mean_sst[b]
}
sst2 <- sst2 %>% dplyr::select(Date, Temperature)


names(sst)[1] <- 'Date'; names(sst2)[1] <- 'Date'; names(wc_sst)[1] <- 'Date'

wc_series <- merge(wc_sst, sst, by='Date')
wc_etuff <- merge(sst, sst2, by='Date')

par(mfrow=c(3,1))
plot(wc_etuff$mean_sst, wc_etuff$Temperature, xlab='WC SST', ylab='eTUFF series')
abline(a=0, b=1)

plot(wc_etuff$min_sst, wc_etuff$Temperature, xlab='WC SST', ylab='eTUFF series')
abline(a=0, b=1)

plot(wc_etuff$max_sst, wc_etuff$Temperature, xlab='WC SST', ylab='eTUFF series')
abline(a=0, b=1)

## CONCLUSION: at least for 1090269, the way we calculate sst from series skews our calculated SSTs much warmer than
## the SSTs found in the -SSTs.csv output from WC which means those SST metrics are calculated differently.

## check another WC tag -> same for 1490108

## check a Lotek tag -> same for 0394



plot(comp_etuff$Temperature, comp_etuff$mean_sst, xlab='eTUFF series SST', ylab='eTUFF sst')
abline(a=0, b=1)
