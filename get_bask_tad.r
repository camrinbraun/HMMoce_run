#dataDir <- '~/ebs/Data/BaskingSharks/batch/'
dataDir <- '~/Documents/WHOI/Data/BaskingSharks/batch/'

meta <- read.table(paste(dataDir, 'bask_metadata.csv',sep=''), sep=',', header=T)
all <- list()

bins <- read.table('~/Documents/WHOI/RData/binQuery.csv', sep=',', header=T)
bins <- bins[which(bins$species == 'BaskingSharks'),]

for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  myDir <- paste(dataDir, ptt, '/', sep='')
  fileList <- list.files(myDir)
  if(length(grep('Histos', fileList)) != 0){
    setwd(myDir)
    tad <- read.table(paste(ptt,'-Histos.csv', sep=''), sep=',', header=T)
    dateVec <- seq(as.Date(meta$TagDate[ii], format='%m/%d/%y'), as.Date(meta$PopDate[ii], format='%m/%d/%y'), by='day')
    
    if (any(tad$HistType == 'TADLimits')){
      # get tadlimits
      
    } else{
      # get them from binquery
      bin.idx <- which(meta[ii,2:4] != '')
      if (bin.idx == 1) tagtype <- 'mk10'
      if (bin.idx == 2) tagtype <- 'mp'
      if (bin.idx == 3) tagtype <- 'mk10af'
      bin.ii <- bins[which(bins$tagyear == meta$TagYear[ii] & bins$tagtype == tagtype),4:ncol(bins)]
    }
    
    tad <- tad[which(tad$HistType == 'TAD'),]
    tad <- tad[,c(2,6,7,14:grep('Bin14', names(tad)))]
    
    tad$Date <- as.Date(tad$Date, format=findDateFormat(tad$Date))
    udates <- unique(tad$Date)
    
    # summarize to one daily set of TAD bins
    for (b in 1:length(udates)){
      tad.b <- tad[which(tad$Date == udates[b]),]
      if (nrow(tad.b) > 1){
        tad.b1 <- tad.b
        tad.b <- tad.b[1,]
        tad.b[,6:ncol(tad.b)] <- colMeans(tad.b1[,6:ncol(tad.b1)])
      }
      
      if (b == 1){
        tad.all <- tad.b
      } else{
        tad.all <- rbind(tad.all, tad.b)
      }
    }
    
    tad <- tad.all
    
    # reshape and assign bin values
    tad <- reshape2::melt(tad, id.vars = c('Ptt','HistType','Date','NumBins','Sum'))
    tad$variable <- tolower(tad$variable)
    tad$bin <- NA
    for (b in 1:length(bin.ii)){
      tad$bin[which(tad$variable %in% names(bin.ii)[b])] <- as.numeric(bin.ii[b])
    }
    
    tad <- tad[which(!is.na(tad$bin)),]
    
    # remove tad data outside of datevec bounds
    tad <- tad[which(tad$Date %in% dateVec),]
    

  }

  all[[ii]] <- tad
  
}
  
  