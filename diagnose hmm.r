meta$hmmoce <- NA
for (i in 37:nrow(meta)){
  setwd(paste('~/ebs/Data/BaskingSharks/batch/', meta$PTT[i],sep=''))
  fileList <- list.files()
  resFiles <- fileList[grep('_res.rda', fileList)]
  if(length(resFiles) == 36){
    meta$hmmoce[i] <- 1
  } else{
    meta$hmmoce[i] <- 0
    print(warning(paste(meta$PTT[i], 'is not complete.')))
  }
  
  if (meta$hmmoce[i] == 1){
    for (ii in 1:length(resFiles)){
      load(resFiles[ii])
      write.table(res$outVec, file='~/ebs/Data/BaskingSharks/batch/baskbatch_res.csv', sep=',', append=T, col.names = F)
    }
  }
  
}

# save to s3
s3save(meta, object='baskbatch_res.rda', bucket='braun-data/Data/BaskingSharks/batch')

all <- read.table('~/ebs/Data/BaskingSharks/batch/baskbatch_res.csv', sep=',', header=F)
names(all) <- c('rownum','ptt','bnd','migr.spd','Lidx','P1','P2','xmin','xmax','ymin','ymax','resolx','resoly','maskL','nll','name')
all <- all[which(all$ptt == ptt),]
all <- all[order(all$nll),]

# once we have the "best" model fit, we can look at diagnostics
fileList <- list.files()
load(fileList[grep('100979_idx7_bndNA_par4', fileList)[1]])
#load(fileList[grep(all$name[nrow(all)], fileList)[1]])
load('check2.rda')

source('~/HMMoce/R/hmm.diagnose.r')
hmm.diagnose(res, plot=TRUE)


