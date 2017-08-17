#=================================
## GET SWORD ENV DATA FOR HMMOCE
#=================================

load('~/Data/Swordfish/batch/unique_sword_dates.rda')
udates <- udatesList; rm(udatesList)
str(udates)
udates <- udates[order(udates)]

# SET SPATIAL LIMITS
sp.lim <- list(lonmin=-85, lonmax=-15,
               latmin=8, latmax=53)

#if (exists('sp.lim')){
#  locs.grid <- setup.locs.grid(sp.lim)
#} else{
#  locs.grid <- setup.locs.grid(locs)
#  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
#                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
#}

# IF YOU NEED TO DOWNLOAD SST DATA
sst.dir <- '~/EnvData/sst/Swordfish/'
dir.create(file.path(sst.dir), recursive = TRUE, showWarnings = FALSE)
udates <- udates[19:length(udates)]
get.env(udates, ptt='sword', type = 'sst', sst.type='ghr', spatLim = sp.lim, save.dir = sst.dir)
get.sst.dates <- udates[!(udates %in% as.Date(substr(list.files(sst.dir), 7, 16)))]
get.env(get.sst.dates, ptt = ptt, type = 'sst', sst.type='ghr',spatLim = sp.lim, save.dir = sst.dir)

# HYCOM DATA
hycom.dir <- '~/EnvData/hycom3/Swordfish/'
dir.create(file.path(hycom.dir), recursive = TRUE, showWarnings = FALSE)
setwd('~/HMMoce'); devtools::load_all()
udates <- udates[order(udates)]
get.env(udates[500], filename='sword', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

get.pdt.dates <- udates[!(udates %in% as.Date(substr(list.files(hycom.dir), 7, 16)))]
get.pdt.dates <- rev(get.pdt.dates)
get.env(get.pdt.dates, filename='sword', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

hycom.dir <- '~/EnvData/hycom3/Swordfish/'
setwd('~/HMMoce'); devtools::load_all()
get.pdt.dates <- udates[!(udates %in% as.Date(substr(list.files(hycom.dir), 7, 16)))]
for (i in 648:1){
  repeat{
    try(get.env(get.pdt.dates[i], filename='sword', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir), silent=T)
    if(file.info(paste('sword_', get.pdt.dates[i],'.nc', sep=''))$size > 37e5){
      break
    }
  }
  print(paste('Finished ', i))
}

save(get.pdt.dates,file='get_pdt_dates.rda')
