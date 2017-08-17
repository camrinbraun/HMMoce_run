#=================================
## GET BASK ENV DATA FOR HMMOCE
#=================================
myDir <- '~/Data/BaskingSharks/batch/'
load('~/Data/BaskingSharks/batch/bask_lims.rda')
str(bask.lims)
#require(HMMoce)
setwd('~/HMMoce')
devtools::load_all()
#----------
# START WITH BIG

# SET SPATIAL LIMITS
sp.lim <- list(lonmin=bask.lims$big.bnds[1], lonmax=bask.lims$big.bnds[2],
               latmin=bask.lims$big.bnds[3], latmax=bask.lims$big.bnds[4])

# IF YOU NEED TO DOWNLOAD SST DATA
sst.dir <- '~/EnvData/sst/BaskingSharks/'
dir.create(file.path(sst.dir), recursive = TRUE, showWarnings = FALSE)
#udates <- udates[19:length(udates)]
get.env(bask.lims$big.dts.sst, filename='bask', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
udates <- bask.lims$big.dts.sst
fileList <- substr(list.files(sst.dir), 10, 19)
fileList <- as.Date(fileList)
get.sst.dates <- udates[!(udates %in% fileList)]

get.sst.dates <- udates[!(udates %in% as.Date(substr(list.files(sst.dir), 6, 15)))]
#get.env(get.sst.dates, ptt = ptt, type = 'sst', sst.type='ghr',spatLim = sp.lim, save.dir = sst.dir)

#setwd()
fileList <- list.files()
for (i in 2:length(fileList)){
  file.rename(fileList[i], paste(substr(fileList[i], 1, 4),'_big', substr(fileList[i], 5, 18),sep=''))
}

# HYCOM DATA
hycom.dir <- '~/EnvData/hycom3/BaskingSharks/'
dir.create(file.path(hycom.dir), recursive = TRUE, showWarnings = FALSE)
udates <- bask.lims$big.dts.pdt
udates <- udates[order(udates)]
fileList <- substr(list.files(hycom.dir), 10, 19)
fileList <- as.Date(fileList[-grep('.nc', fileList)])
get.pdt.dates <- udates[!(udates %in% fileList)]

get.env(udates, filename='bask', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
fileList <- list.files()

#for (i in 2:length(fileList)){
#  file.rename(fileList[i], paste(substr(fileList[i], 1, 4),'_big', substr(fileList[i], 5, 18),sep=''))
#}

#===============
# NOW DO SMALLER SPATIAL EXTENT

# SET SPATIAL LIMITS
sp.lim <- list(lonmin=bask.lims$small.bnds[1], lonmax=bask.lims$small.bnds[2],
               latmin=bask.lims$small.bnds[3], latmax=bask.lims$small.bnds[4])

# IF YOU NEED TO DOWNLOAD SST DATA
sst.dir <- '~/EnvData/sst/BaskingSharks/'
#dir.create(file.path(sst.dir), recursive = TRUE, showWarnings = FALSE)
#udates <- udates[19:length(udates)]
get.env(bask.lims$small.dts.sst, filename='bask_small', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
udates <- bask.lims$small.dts.sst
get.sst.dates <- udates[!(udates %in% as.Date(substr(list.files(sst.dir), 12, 21)))]
#get.env(get.sst.dates, ptt = ptt, type = 'sst', sst.type='ghr',spatLim = sp.lim, save.dir = sst.dir)

#setwd()
fileList <- list.files()
for (i in 2:length(fileList)){
  file.rename(fileList[i], paste(substr(fileList[i], 1, 4),'_small', substr(fileList[i], 5, 18),sep=''))
}

# HYCOM DATA
hycom.dir <- '~/EnvData/hycom3/BaskingSharks/'
#dir.create(file.path(hycom.dir), recursive = TRUE, showWarnings = FALSE)
udates <- bask.lims$small.dts.pdt
udates <- udates[order(udates)]
get.pdt.dates <- udates[!(udates %in% as.Date(substr(list.files(hycom.dir), 6, 15)))]
get.env(get.pdt.dates[2:length(get.pdt.dates)], filename='bask', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
fileList <- list.files()

#for (i in 2:length(fileList)){
#  file.rename(fileList[i], paste(substr(fileList[i], 1, 4),'_small', substr(fileList[i], 5, 18),sep=''))
#}
