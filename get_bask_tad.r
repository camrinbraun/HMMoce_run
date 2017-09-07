dataDir <- '~/ebs/Data/BaskingSharks/batch/'
#dataDir <- '~/Documents/WHOI/Data/BaskingSharks/batch/'

meta <- read.table(paste(dataDir, 'bask_metadata.csv',sep=''), sep=',', header=T)
meta <- meta[which(!is.na(meta$res)),]

all <- list()

bins <- read.table(paste(dataDir, 'binQuery.csv', sep=''), sep=',', header=T)
bins <- bins[which(bins$species == 'BaskingSharks'),]

bin.c.all <- list(list(1,2,3,4:7,8:9,10:12, NA),
                  list(1,2,3,4:7,8:9,10:12, NA),
                  list(1,2,3,4:7,8:9,10:12, NA),
                  list(1:2,3,4,5:6,7:8,9:13,14),
                  list(1,2,3,4:7,8,9:11,12),
                  list(1:2,3,4,5:7,8,9:11,12),
                  list(1:2,3,4,5:8,9:10,11:13,14),
                  list(1:2,3,4,5:8,9:10,11:13,14))
#lapply(bin.c.all, FUN=function(x) length(x))

for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  myDir <- paste(dataDir, ptt, '/', sep='')
  fileList <- list.files(myDir)
  if(length(grep('Histos', fileList)) != 0){
    setwd(myDir)
    tad <- read.table(paste(ptt,'-Histos.csv', sep=''), sep=',', header=T)
    dateVec <- seq(as.Date(meta$TagDate[ii], format='%m/%d/%y'), as.Date(meta$PopDate[ii], format='%m/%d/%y'), by='day')
    
    #if (any(tad$HistType == 'TADLimits')){
      # get tadlimits
      
    #} else{
      # get them from binquery
      bin.idx <- which(meta[ii,2:4] != '')
      if (bin.idx == 1) tagtype <- 'mk10'
      if (bin.idx == 2) tagtype <- 'mp'
      if (bin.idx == 3) tagtype <- 'mk10af'
      bin.ii <- bins[which(bins$tagyear == meta$TagYear[ii] & bins$tagtype == tagtype),4:ncol(bins)]
      bin.c <- bin.c.all[[which(bins$tagyear == meta$TagYear[ii] & bins$tagtype == tagtype)]]
    #}
    
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
      
      sel <- tad.b[,1:12]
      for (r in 1:length(bin.c)){
        if (is.na(bin.c[[r]])){
          sel[,r+5] <- 0
        } else{
          sel[,r+5] <- sum(tad.b[,c(bin.c[[r]]+5)], na.rm=T)
        }
      }
      sel$Sum <- sum(sel[,6:12], na.rm=T)
      
      if (b == 1){
        tad.all <- sel
      } else{
        tad.all <- rbind(tad.all, sel)
      }
    }
    
    tad <- tad.all
    
    # reshape and assign bin values
    tad <- reshape2::melt(tad, id.vars = c('Ptt','HistType','Date','NumBins','Sum'))
    tad$variable <- tolower(tad$variable)
    tad$bin <- NA
    bins.all <- data.frame(matrix(c(10,25,50,200,400,1000,2000), ncol=7, byrow=T))
    names(bins.all) <- names(bin.ii)[1:7]
    for (b in 1:length(bins.all)){
      tad$bin[which(tad$variable %in% names(bins.all)[b])] <- as.numeric(bins.all[b])
    }
    
    tad <- tad[which(!is.na(tad$bin)),]
    
    # remove tad data outside of datevec bounds
    tad <- tad[which(tad$Date %in% dateVec),]
    
    #tad <- tad[which(tad$value != 0),]

  }

  all[[ii]] <- tad
  rm(bin.c)
}
  
for (b in 1:length(all)){
  if (b == 1){
    df <- all[[b]]
  } else{
    df <- rbind(df, all[[b]])
  }
  
}
df <- df[,c(1,3,6:8)]
names(df) <- tolower(names(df))

df$season <- NA
for (b in 2004:2012){
  # WINTER
  if(b != 2004) df$season[which(df$date >= paste(b-1,'-12-20',sep='') & df$date <= paste(b-1,'-12-31',sep=''))] <- 4
  df$season[which(df$date >= paste(b,'-01-01',sep='') & df$date < paste(b,'-03-20',sep=''))] <- 4
  # FALL
  df$season[which(df$date >= paste(b,'-09-20',sep='') & df$date < paste(b,'-12-20',sep=''))] <- 3
  #SPR
  df$season[which(df$date >= paste(b,'-03-20',sep='') & df$date < paste(b,'-06-20',sep=''))] <- 1
  #SUM
  df$season[which(df$date >= paste(b,'-06-20',sep='') & df$date < paste(b,'-09-20',sep=''))] <- 2
}

#library(dplyr)
summ <- data.frame(group_by(df, variable, season) %>% summarise(mean(value, na.rm=T)))

setwd(dataDir)
all.tad.res <- list(alltadList=all, bigdf=df, summ=summ)
save(all.tad.res, file='all_bask_tad_v3.rda')

  