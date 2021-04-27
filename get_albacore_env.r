devtools::install_github('camrinbraun/HMMoce', ref='dev')
library(HMMoce)
library(tidyverse)
source('../analyzePSAT/R/make360.R')
library(fields) ## for quick mapping
library(raster)

base_dir <- '~/ebs/RCode/HMMoce_run'
setwd(base_dir)

meta <- read.table('../nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, stringsAsFactors = F)
meta <- meta %>% filter(platform == 'Thunnus alalunga' & etuff == 1)
meta$time_coverage_start <- as.POSIXct(meta$time_coverage_start, tz = 'UTC')
meta$time_coverage_end <- as.POSIXct(meta$time_coverage_end, tz = 'UTC')

## a set of manually defined spatial bounds for each fish
limits <- read.table('~/ebs/Data/albacore/alb_limits.csv', sep=',', header=T)

#for (i in 1:9){
for (i in 1:nrow(meta)){
  
  
  ## temporal bounds
  iniloc <- data.frame(matrix(c(lubridate::day(meta$time_coverage_start[i]),
                                lubridate::month(meta$time_coverage_start[i]),
                                lubridate::year(meta$time_coverage_start[i]),
                                meta$geospatial_lat_start[i],
                                meta$geospatial_lon_start[i],
                                lubridate::day(meta$time_coverage_end[i]),
                                lubridate::month(meta$time_coverage_end[i]),
                                lubridate::year(meta$time_coverage_end[i]),
                                meta$geospatial_lat_end[i],
                                meta$geospatial_lon_end[i]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
  tag <- iniloc$date[1]
  pop <- iniloc$date[2]
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- seq.POSIXt(tag, pop, by = '24 hours')

  sp.lim <- list(lonmin = make360(limits$lonmin[which(limits$instrument_name == meta$instrument_name[i])]),
                 lonmax = make360(limits$lonmax[which(limits$instrument_name == meta$instrument_name[i])]),
                 latmin = limits$latmin[which(limits$instrument_name == meta$instrument_name[i])],
                 latmax = limits$latmax[which(limits$instrument_name == meta$instrument_name[i])])
  
  ## setup the spatial grid to base likelihoods on
  locs.grid <- setup.locs.grid(sp.lim, res='quarter')
  
  work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
  if (!dir.exists(work_dir)) dir.create(work_dir)
  setwd(work_dir)
  
  sst.dir <- paste0('./tmp/sst/')
  if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
  for (b in 1:length(dateVec)){
    if (!file.exists(paste0(sst.dir, '/', 'oi_', dateVec[b], '.nc'))) get.env(dateVec[b], filename='oi', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
  }
  
  hycom.dir <- paste0('./tmp/hycom/')
  if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
  
  ## check for hycom. if its not available pull it down from the s3 bucket
  #if (length(list.files(hycom.dir)) != length(dateVec)){
  #  system(paste0('aws s3 cp s3://braun-data/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/ /home/rstudio/ebs/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/'))
  #  if (length(list.files(hycom.dir)) != length(dateVec)) stop('Check that all daily hycom data is available for this individual run.')
  #}
  
  ## this is how we get hycom. it will skip days for which data already exists
  for (b in 1:length(dateVec)){
    setwd(work_dir)
    if (!file.exists(paste0(hycom.dir, '/', 'hycom_', dateVec[b], '.nc'))) {
      try(
        get.env(dateVec[b], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir), TRUE)
      
      #if (b %% 10 == 0){
      #removeTmpFiles(); 
      removeTmpFiles(h=0)
      file.remove(list.files(tempdir(), recursive=TRUE))
      #gc()
      #}
    }
  }
  
}


i=9
paste0('aws s3 sync ./', meta$instrument_name[i], '/ s3://braunmpg/RCode/HMMoce_run/', meta$instrument_name[i], '/')# --exclude "*" --include "*.nc" --dryrun')


for (i in 2:9){
  fList <- list.files(paste0('~/ebs/RCode/HMMoce_run/', meta$instrument_name[i], '/tmp/hycom/'), recursive=TRUE, full.names = TRUE)
  file.remove(fList)
  
}
