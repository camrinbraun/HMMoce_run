
dataDir <- '~/ebs/Data/BaskingSharks/batch/'
meta <- read.table(paste(dataDir, 'bask_metadata.csv',sep=''), sep=',', header=T)

for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  setwd(paste('~/ebs_small/Data/BaskingSharks/batch/', ptt, '/', sep='')) 
  homeDir <- paste('~/ebs/Data/BaskingSharks/batch/', ptt, '/', sep='')
  
  fileList <- list.files()
  file.idx <- grep('check2.rda', fileList)
  if (length(file.idx) == 1){
    file.copy(fileList[file.idx], paste(homeDir, fileList[file.idx],sep=''))
    
    #get.pdfs <- grep('_results.pdf', fileList)
    #file.copy(fileList[file.idx], paste(homeDir, fileList[file.idx],sep=''))
    #file.copy(fileList[get.pdfs], paste(homeDir, fileList[get.pdfs],sep=''))
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
