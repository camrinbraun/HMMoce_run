#setwd('~/Documents/WHOI/RCode/HMMoce/'); devtools::load_all()
dataDir <- '~/ebs/Data/BaskingSharks/batch/'

meta <- read.table(paste(dataDir, 'bask_metadata.csv',sep=''), sep=',', header=T)
#meta <- meta[-5,]
all <- list()
for (ii in 1:nrow(meta)){
  ptt <- meta$PTT[ii]
  myDir <- paste(dataDir, ptt, '/', sep='')
  # load individual results, compile into a list
  if(!is.na(meta$res[ii])){
    load(paste(myDir, meta$res[ii],'-HMMoce_res.rda', sep=''))
    
    # run getCtr
    bnds <- getCtr(res$s, res$tr, res$g, threshold = 10, makePlot=F)
    
    # compile df of position, date, bnds (xtracto)
    #data.frame(lapply(bnds, FUN=function(x) c(x$yDist, x$xDist)))
    df <- cbind(res$tr, t(data.frame(lapply(bnds, FUN=function(x) c(x$yDist, x$xDist)))))
    names(df)[5:6] <- c('ydist','xdist')
    df$ydist[which(is.na(df$ydist))] <- mean(df$ydist, na.rm=T)
    df$xdist[which(is.na(df$xdist))] <- mean(df$xdist, na.rm=T)
    df$ptt <- ptt
    rownames(df) <- NULL
    
    # make list of each RD
    rd <- plotRD(res$s, res$tr, ptt, g=res$g, makePlot=F)
    all[[ii]] <- list(allRD=rd$allRD, behavRD=rd$behavRD, df=df)
    
  }
  
}

# remove NULL elements
names(all) <- meta$PTT
all[sapply(all, is.null)] <- NULL

#lapply(all, FUN=function(x) which(x$allRD))
allRD <- list()
rsamp <- all[[5]]$allRD
rsamp <- disaggregate(rsamp, 4)

for (b in 1:length(all)){
  if(!is.null(all[[b]])){
    allRD[[b]] <- all[[b]]$allRD
    allRD[[b]] <- resample(allRD[[b]], rsamp)
  }
  
}


#names(allRD) <- meta$PTT
#allRD[sapply(allRD, is.null)] <- NULL


s <- stack(allRD)
#pdf('sword_map.pdf', width=12, height=8)
plot(ln(sum(s, na.rm=T)))
world(add=T, fill=T)

# get total df
rb <- lapply(all, FUN=function(x) x$df)
for (i in 1:length(rb)){
  if(i==1){
    all.df <- rb[[i]]
  } else{
    all.df <- rbind(all.df, rb[[i]])
  }
  lines(rb[[i]]$lon, rb[[i]]$lat)
}

bask.res <- list(all=all, allRD=s, all.df=all.df)
save(bask.res, file='bask_results.rda')

#dev.off()
#allRDs <- sum(allRD, na.rm=T)#lapply(allRD, FUN=function(x) sum(x, na.rm=T))

#allRD <- lapply(all, FUN=function(x) brick(x$allRD))

#=========================
## END


# skip 2, 3 for now
pttList <- c(98721, 100980, 104668, 104671, 104672, 106788, 110490, 110491)
for (ii in 1:length(pttList)){
  ptt <- pttList[ii]
  dataDir <- '~/Documents/WHOI/Data/Swordfish/batch/'
  myDir <- paste(dataDir, ptt, '/', sep='')
  # load individual results, compile into a list
  load(paste(myDir, ptt,'-HMMoce.RData', sep=''))
  dataDir <- '~/Documents/WHOI/Data/Swordfish/batch/'
  str(g)
  res <- list(outVec = outVec, s = s, tr = tr, g=g, dateVec=dateVec, iniloc=iniloc, grid = raster::res(L.rasters[[resamp.idx]])[1])
  str(res)
  #invisible(readline(prompt="Press [enter] to perform the next iteration and plot"))
  base::save(res, file=paste(dataDir, ptt, '/', ptt, '-HMMoce_res1.RData', sep=''))
  rm(list = ls())
  
  
}
