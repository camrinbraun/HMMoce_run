
# SETWD
setwd('~/Documents/WHOI/RCode/HMMoce_run/data/141259/') 

# READ IN TAG DATA
ptt <- 141259

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 
                              10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
#myDir <- '~/Documents/WHOI/RCode/HMMoce/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
myDir <- getwd()
# sst data
tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

# light data
#light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop);
light <- read.table('141259-LightLoc.csv', sep=',', header=T, skip=2)
light.udates <- unique(as.Date(light$Day, format='%d-%b-%y'))#; light <- light$data
light.udates <- light.udates[!is.na(light.udates)]
light.udates <- light.udates[which(light.udates >= as.Date(tag) & light.udates <= as.Date(pop))]
light.udates <- light.udates[3:length(light.udates)]

#blue254.light <- sample(light.udates, round(.25 * length(light.udates), 0))
#blue254.light <- blue254.light[order(blue254.light)]
#blue254.sst <- sample(sst.udates, round(.5 * length(sst.udates), 0))
#blue254.sst <- blue254.sst[order(blue254.sst)]

#blue256.light <- sample(light.udates, round(.25 * length(light.udates), 0))
#blue256.light <- blue256.light[order(blue256.light)]
#blue256.sst <- sample(sst.udates, round(.5 * length(sst.udates), 0))
#blue256.sst <- blue256.sst[order(blue256.sst)]

blue259.light <- sample(light.udates, round(.25 * length(light.udates), 0))
blue259.light <- blue259.light[order(blue259.light)]
blue259.sst <- sample(sst.udates, round(.5 * length(sst.udates), 0))
blue259.sst <- blue259.sst[order(blue259.sst)]
