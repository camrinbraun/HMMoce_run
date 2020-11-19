## script to clean up the albacore tags

home_dir <- '~/Google Drive File Stream/My Drive/Albacore - All Data/Tags w o Migrations/'
dirList <- list.files(home_dir)
dirList <- dirList[which(dirList != 'Icon\r')]

meta <- read.table('~/Google Drive File Stream/My Drive/Albacore - All Data/Albacore Tag Summary Mar2020 CB_etuff.csv', sep=',', header=T, blank.lines.skip = F, stringsAsFactors = F)
#meta <- meta %>% filter(!is.na(Days.free))
#meta <- meta[,1:which(names(meta) %in% 'Days.free')]
meta$time_coverage_end <- as.Date(meta$time_coverage_end, format='%m/%d/%y')
meta$time_coverage_start <- as.Date(meta$time_coverage_start, format='%m/%d/%y')

meta$geospatial_lon_end <- make360(meta$geospatial_lon_end)
meta$geospatial_lon_start <- make360(meta$geospatial_lon_start)


####################### "A2082 - Lotek"  #######################

id <- 'A2082'
idx <- which(meta$serial_number %in% id)

data_dir <- paste(home_dir, dirList[which(dirList == id)], '/', sep='')
fList <- list.files(data_dir, full.names = T)
data_file <- fList[grep('TIMESERIES', fList)]
summary_file <- fList[grep('DAYLOG', fList)]
df <- data.table::fread(data_file, sep=',', header=T)#, nrows=Inf)
df$dt <- as.POSIXct(df$`Time (UTC)`, format='%m/%d/%Y %H:%M:%S', tz='UTC')

start <- as.POSIXct(meta$Release.date[idx], tz='UTC')
end <- as.POSIXct(meta$End.date[idx], tz='UTC')

p1 <- ggplot(df, aes(x=dt, y=`Depth - dBar`, colour=`Ext Temp deg C`)) +
  geom_point() +
  geom_vline(xintercept = start) +
  geom_vline(xintercept = end)


png(paste('~/Google Drive File Stream/My Drive/Albacore - All Data/figures/', id, '-depthtemp_raw.png', sep=''), width=12, height=8, units='in', res=300)
p1
dev.off()

## filter
df <- df %>% filter(dt <= end)

p1 <- ggplot(df, aes(x=dt, y=`Depth - dBar`, colour=`Ext Temp deg C`)) +
  geom_point() +
  geom_vline(xintercept = start) +
  geom_vline(xintercept = end)
png(paste('~/Google Drive File Stream/My Drive/Albacore - All Data/figures/', id, '-depthtemp_clean.png', sep=''), width=12, height=8, units='in', res=300)
p1
dev.off()


p1 <- ggplot(df, aes(x=dt, y=`Ext Temp deg C`)) +
  geom_line(colour='blue') +
  geom_line(data = df, aes(x=dt, y=`Int Temp deg C`), colour='red') #+
  #geom_vline(xintercept = as.POSIXct('2005-05-10')) #+
  #annotate('text', label = 'external thermistor goes bad', x=as.POSIXct('2005-05-10'), y=50)
png(paste('~/Google Drive File Stream/My Drive/Albacore - All Data/figures/', id, '-comparetemp.png', sep=''), width=12, height=8, units='in', res=300)
p1
dev.off()

### Clean summary

summ <- data.table::fread(summary_file, sep=',', header=T, nrows=Inf)
summ$dt <- as.POSIXct(summ$`Mission Date`, format='%m/%d/%Y', tz='UTC')
summ <- summ %>% filter(dt >= start & dt <= end)
names(summ) <- c('date','sunrise_utc','sunset_utc','lon_txt','lat_txt','sst_min','sst_max','sst_median', 'dt')

summ$lon <- NA
summ$lon <- as.numeric(summ$lon)
summ$lat <- NA
summ$lat <- as.numeric(summ$lat)

for (i in 1:nrow(summ)){
  ## LONGITUDE
  summ$lon[i] <- as.numeric(substr(summ$lon_txt[i], 
                                1, 
                                stringr::str_locate_all(summ$lon_txt[i], ' ')[[1]][1]))
  is_west <- ifelse(substr(summ$lon_txt[i], 
                           stringr::str_locate_all(summ$lon_txt[i], ' ')[[1]][1] + 1,
                           nchar(summ$lon_txt)) == 'West',
                    TRUE, FALSE)
  if (is_west) summ$lon[i] <- summ$lon[i] * -1
  
  
  ## LATITUDE
  summ$lat[i] <- as.numeric(substr(summ$lat_txt[i], 
                                   1, 
                                   stringr::str_locate_all(summ$lat_txt[i], ' ')[[1]][1]))
  is_south <- ifelse(substr(summ$lat_txt[i], 
                           stringr::str_locate_all(summ$lat_txt[i], ' ')[[1]][1] + 1,
                           nchar(summ$lat_txt)) == 'South',
                    TRUE, FALSE)
  if (is_south) summ$lat[i] <- summ$lat[i] * -1
                         
}



iniloc <- data.frame(matrix(c(lubridate::day(start), lubridate::month(start), lubridate::year(start), meta$geospatial_lon_start[idx], meta$geospatial_lat_start[idx],
                              lubridate::day(end), lubridate::month(end), lubridate::year(end), meta$geospatial_lon_end[idx], meta$geospatial_lat_end[idx]), nrow = 2, ncol = 5, byrow = T))
names(iniloc) <- list('day','month','year','lat','lon')
iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
tag <- iniloc$date[1]
pop <- iniloc$date[2]

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- seq.POSIXt(tag, pop, by = '24 hours')


## MMD
## create max depth df for calc.bathy
df$Date <- as.Date(df$dt, tz='UTC')
mmd <- data.frame(df %>% group_by(Date) %>% summarise(MaxDepth = max(`Depth - dBar`)))
mmd$Date <- as.POSIXct(mmd$Date, tz='UTC')
mmd$MaxDepth <- mmd$MaxDepth * -1
mmd <- mmd[which(mmd$Date <= iniloc$date[2]),]
udates <- unique(as.Date(mmd$Date))
gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have MMD data...', sep=''))


## PDT
## create depth-temperature profile data for various 3d likelihood calcs
pdt_dates <- seq.POSIXt(tag, pop, by = 'day') ## generate daily summary of depth-temp profiles
pdt <- df[,c('dt','Depth - dBar','Ext Temp deg C')]
names(pdt) <- c('Date','Depth','Temperature')
pdt$Depth <- pdt$Depth * -1
## generate depth-temp summary from time series
pdt <- bin_TempTS(pdt, out_dates = pdt_dates, bin_res = 25)
pdt <- pdt[,c('Date','Depth','MinTemp','MaxTemp')]
udates <- unique(as.Date(pdt$Date))
gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have PDT data...', sep=''))


## SST
## assume daily max temp is equivalent to SST
sst <- data.frame(summ[,c('dt','sst_median')])
names(sst) <- c('Date','Temperature')
udates <- unique(as.Date(sst$Date))
gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have SST data...', sep=''))


## LIGHT LOCS
#### speed filter to eliminate spurious locs?
locs <- data.frame(summ)
locs <- locs[,c('dt','lat','lon')]
locs <- locs[order(locs$dt),]
names(locs) = c('Date','Latitude','Longitude')
locs$Error.Semi.minor.axis = .7 * 1000 * 111

{
  time_step <- 24
  res_out <- 24
  
  summ_df <- summ
  summ_df$id <- id
  summ_df$lc <- 'Z'
  summ_df <- summ_df[,c('id','dt','lc','lon','lat')]
  names(summ_df)[2] <- 'date'
  summ_df[nrow(summ_df) + 1,] <- c(id, start, 'Z', meta$geospatial_lon_start[idx], meta$geospatial_lat_start[idx])
  summ_df[nrow(summ_df) + 1,] <- c(id, end, 'Z', meta$geospatial_lon_end[idx], meta$geospatial_lat_end[idx])
  summ_df <- summ_df[order(summ_df$date),]
  
  df.locs <- split(summ_df, summ_df$id)
  
  ssm_fit <- lapply(df.locs, function(x) foieGras::fit_ssm(x, model='crw', time.step = time_step, vmax=10))
  
  # which did or did not converge
  idx <- lapply(ssm_fit, function(x) x$converged) %>% do.call(rbind, .)
  
  if (any(!idx)){
    for (ii in which(!idx)){
      print('Some initial use of fit_ssm() did not converge. Trying other time_step values.')
      #df.ii <- foieGras::grab(ssm_fit[[ii]], 'data', as_sf=F)
      df.ii <- df.locs[[ii]]
      #df.ii <- ssm_fit[[ii]]$ssm[[1]]$data
      df.ii <- data.frame(df.ii[,names(df)])
      
      # iterate through some reasonable time steps
      poss.tstep <- c(6,3,1,12,24)
      poss.tstep <- poss.tstep[-which(poss.tstep %in% time_step)]
      for (bb in poss.tstep){
        try_fit <- foieGras::fit_ssm(df.ii, model='crw', time.step = bb, vmax=10, optim='nlminb')
        if (try_fit$converged){
          ssm_fit[[ii]] <- try_fit
          print(paste('Converged on time step', bb, 'hours.'))
          break
        } else if (!try_fit$converged & bb == 24){
          print('None of the alternative time steps resulted in fit_ssm() model convergence.')
          print('Cancelling simulations.')
          sim <- FALSE; tmp_out <- NULL
        } else{
          next
        }
      } # bb
    } # which idx
  } # if any
  
  #ssm_fit <- foieGras::fit_ssm(df, model='crw', time_step = time_step, vmax=10)
  t2 <- Sys.time()
  
  idx <- lapply(ssm_fit, function(x) x$converged) %>% do.call(rbind, .)
  ## drop those trajectories that still don't converge
  if (any(!idx)){
    ssm_fit <- ssm_fit[which(idx)]
    
    ## if dropping trajs that dont converge result in no data left, stop this run
    if (length(ssm_fit) == 0) stop('fit_ssm did not converge despite trying several time steps between 1 and 24 hours. Check the input data for errors.')
    
    cat(sprintf("dropping %d trajectories due to lack of convergence in fit_ssm(). %d trajectories remain.\n", length(which(!idx)), length(which(idx))))
    
  }
  
  cat(sprintf("\nfit took %d seconds...\n", round(as.numeric(difftime(t2, t1, units = 'secs')), 0)))
  
  ## grab predicted locations output from fit_ssm
  plocs <- lapply(ssm_fit, FUN=function(x) foieGras::grab(x, what = "p", as_sf = FALSE)) %>%
    do.call(rbind, .) %>%
    tbl_df() %>%
    mutate(id = as.character(id)) %>% group_by(id)
  
  ## subsample from predicted locations, if applicable
  if (!is.null(res_out)){
    plocs <- split(plocs, plocs$id)
    
    for (i in 1:length(plocs)){
      res.i <- as.numeric(difftime(plocs[[i]]$date[2], plocs[[i]]$date[1], units = 'hours'))
      res_slice <- res_out / res.i
      if (res_slice < 1) res_slice <- 1
      plocs[[i]] <- plocs[[i]] %>% do(slice(., seq(1, n(), by = res_slice)))
    }
    
    plocs <- plocs %>%
      do.call(rbind, .) %>%
      tbl_df() %>%
      mutate(id = as.character(id)) %>% group_by(id)
    
  }
  
}


sp.lim <- list(lonmin = -145,
               lonmax = -110,
               latmin = 15,
               latmax = 40)

## setup the spatial grid to base likelihoods on
locs.grid <- setup.locs.grid(sp.lim, res='quarter')

udates <- seq.Date(as.Date(tag), as.Date(pop), by = 'day')

dir <- paste0('~/work/RCode/HMMoce_run/', id, '/')
if (!dir.exists(dir)) dir.create(dir)
setwd(dir)
#load('./141259_tryHMMoce_20200724.rda')
#load('./141259_L.res_20200727.rda')

sst.dir <- paste0(dir, '/EnvData/sst/')
if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
for (i in 1:length(udates)) get.env(udates[i], filename='oisst', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)

hycom.dir <- paste0(dir,'/EnvData/hycom/')
if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
for (i in 1:length(udates)) get.env(udates[i], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

bathy.dir <- paste0(dir, '/EnvData/bathy/')
if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
bathy <- HMMoce::get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir, res=1)
} else{ ## OR (once downloaded and reading the .nc later)
bathy <- irregular_ncToRaster(paste0(bathy.dir, 'bathy.nc'), varid = 'topo')
}


L.light <- calc.lightloc(locs, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE)

L.sst <- calc.sst(sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1, focalDim = 3)


# OCEAN HEAT CONTENT (INTEGRATED PDTS)
L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
#L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)

# HYCOM PROFILE BASED LIKELIHOODS
L.hycom.se <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
#L.hycom.par <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:5], use.se = F, ncores=2)

bathy_resamp <- raster::resample(bathy, L.sst) # or whatever grid makes sense to resample to
L.bathy <- calc.bathy(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')

## make list of rasters
L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.hycom.se = L.hycom.se, L.bathy = L.bathy)

## resample rasters
resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
save(L.res, file=paste0(id, '_L.res_', Sys.Date(), '.rda'))
save(L.rasters, paste0(id, '_L.rasters_', Sys.Date(), '.rda'))










