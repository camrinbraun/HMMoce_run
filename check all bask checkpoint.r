
dataDir <- '~/ebs/Data/BaskingSharks/batch/'
meta <- read.table(paste(dataDir, 'bask_metadata.csv',sep=''), sep=',', header=T)

for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  setwd(paste('~/ebs_small/Data/BaskingSharks/batch/', ptt, '/', sep='')) 
  homeDir <- paste('~/ebs/Data/BaskingSharks/batch/', ptt, '/', sep='')
  
  fileList <- list.files()
  file.idx <- grep('_res.rda', fileList)
  if (length(file.idx) == 36){
    file.copy(fileList[file.idx], paste(homeDir, fileList[file.idx],sep=''))
    
    get.pdfs <- grep('_results.pdf', fileList)
    #file.copy(fileList[file.idx], paste(homeDir, fileList[file.idx],sep=''))
    file.copy(fileList[get.pdfs], paste(homeDir, fileList[get.pdfs],sep=''))
    print(paste('Files copied for ', ptt, '...', sep=''))
    
  } else if (length(file.idx) > 36){
    warning(paste('Error: file index longer than expected. Please check ', ptt, '...', sep=''))
  } else if (length(file.idx) == 0){
    #print(paste('No files exist for ', ptt, '...', sep=''))
    #print(paste('Trying home directory...'))
    fileList <- list.files(homeDir)
    file.idx <- grep('check2.rda', fileList)
    if (length(file.idx) == 1){
      print(paste('Files already home for ', ptt, '...', sep=''))
    } else {
      warning(paste('Error: No files home for ', ptt, '...', sep=''))
    }
  } else{
    warning(paste('Something strange with ', ptt, '...'))
  }
  
  
}

#ptt <- meta$PTT[ii]
setwd(paste('~/ebs_small/Data/BaskingSharks/batch/', ptt, '/', sep='')) 
list.files()
homeDir <- paste('~/ebs/Data/BaskingSharks/batch/', ptt, '/', sep='')
list.files()

for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  setwd(paste('~/ebs_small/Data/BaskingSharks/batch/', ptt, '/', sep='')) 
  #homeDir <- paste('~/ebs/Data/BaskingSharks/batch/', ptt, '/', sep='')
  
  fileList <- list.files()
  file.idx <- grep('_res.rda', fileList)
  file.remove(fileList[file.idx])
  file.idx <- grep('_results.pdf', fileList)
  file.remove(fileList[file.idx])
  
  
}


  if(length(file.idx) == 1){
    print(paste(ptt, ' complete'))
  } else{
    print(paste(ptt, ' NOT'))
  }
}






dataDir <- '~/ebs/Data/BaskingSharks/batch/'
meta <- read.table(paste(dataDir, 'bask_metadata.csv',sep=''), sep=',', header=T)
naList <- c(67812, 88144, 88145, 95979, 95978)
meta <- meta[-which(meta$PTT %in% naList),]

for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  setwd(paste('~/ebs/Data/BaskingSharks/batch/', ptt, '/', sep='')) 
  
  fileList <- list.files()
  file.idx <- grep('_res.rda', fileList)
  
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
  #load(paste(myDir, ptt,'_likelihoods2.RData', sep=''))
  
  # sst data
  tag.sst <- read.wc2(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
  #sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data
  
  # depth-temp profile data
  require(dplyr)
  pdt <- read.wc2(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
  #pdt.udates <- pdt$udates; pdt <- pdt$data
  
  # use GPE2 locs for light-based likelihoods
  locs <- read.table(paste(myDir, ptt, '-Locations-GPE2.csv', sep = ''), sep = ',', header = T, blank.lines.skip = F)
  locDates <- as.Date(as.POSIXct(locs$Date, format=HMMoce:::findDateFormat(locs$Date)))
  gpe.count <- length(which(as.Date(seq(tag, pop, 'day')) %in% unique(locDates)))
  gpe.gaps <- diff(c(as.Date(tag), unique(locDates), as.Date(pop)), units='days')
  gpe.gaps <- gpe.gaps[gpe.gaps > 1]
  
  for (b in file.idx){
    load(fileList[b])
    outVec.all <- c(unlist(res$outVec), tag.sst$deploydur, tag.sst$data.count,
                    min(tag.sst$gaps), max(tag.sst$gaps), mean(tag.sst$gaps),
                    pdt$data.count, min(pdt$gaps), max(pdt$gaps), mean(pdt$gaps),
                    gpe.count, min(gpe.gaps), max(gpe.gaps), mean(gpe.gaps))
    write.table(outVec.all, file='~/ebs/Data/BaskingSharks/batch/all_bask_outvec.csv', sep=',', col.names=F, row.names=T, append=T)
    
  }
  
}


all1 <- read.table('~/ebs/Data/BaskingSharks/batch/all_bask_outvec.csv', sep=',', header=F)

for (i in 17:nrow(meta)){
  ptt <- meta$PTT[i]
  idx <- which(all[,2] == ptt)
  
  for (b in 1:length(idx)){
    getall <- data.frame(t(as.character(all[idx[b]:(idx[b]+27),2])))
    write.table(getall, '~/ebs/Data/BaskingSharks/batch/all_bask_outvec.csv', sep=',', col.names=F, append=T)
  }
}

  if (i == 1){
    allt <- data.frame(matrix(t(all[1:28,2]), ncol=28))
  } else{
    idx <- i:(i+27)
    allt <- rbind(allt, data.frame(matrix(t(all[1:28,2]), ncol=28)))
  }
  
  
}
