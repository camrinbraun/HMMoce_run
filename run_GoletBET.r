devtools::load_all('~/work/RCode/HMMoce')
library(raster)
#library(foieGras)
library(tidyverse)

#meta <- read.table('../nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, stringsAsFactors = F)
xlsFiles <- list.files('~/work/Data/Bigeye/WaltTagsAug2020/', full.names = T, recursive = T)
xlsFiles <- xlsFiles[grep('.xls', xlsFiles)]
xlsFiles <- xlsFiles[1]

#for (aa in 1:length(xlsFiles)){
  aa = 1
  xlsfile = xlsFiles[aa]
  locs <- gdata::read.xls(xlsfile, sheet='Lat&Long', skip=1, header=T, stringsAsFactors=F)#[,1:6] # 6 cols
  # dat = MWTextract(tagID = 57508, xlsfile, delta = T, minmax = F)
  day0 = as.POSIXct(locs[3,8], format='%b %d, %Y', tz='UTC') ## tag date
  x0 = as.numeric(locs[5,8:9])#(read_excel(xlsfile, skip = 5, sheet = 'Lat&Long', n_max = 1))[,8:9])
  x0[2] <- x0[2] * -1
  dayT = as.POSIXct(locs[9,8], format='%b %d, %Y', tz='UTC') ## end date
  xT = as.numeric(locs[11,8:9])
  xT[2] <- xT[2] * -1
  
  iniloc <- data.frame(matrix(c(lubridate::day(day0), lubridate::month(day0), lubridate::year(day0), x0[1], x0[2],
                                lubridate::day(dayT), lubridate::month(dayT), lubridate::year(dayT), xT[1], xT[2]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
  tag <- iniloc$date[1]
  pop <- iniloc$date[2]
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  
  ## read MT tag data
  light <- gdata::read.xls(xlsfile, sheet='Sunrise and Sunset Times', skip=1, header=T)[,1:5]
  #lightloc <- gdata::read.xls(xlsfile, sheet='Lat&Long', skip=1, header=T)#[,1:6] # 6 cols
  depth <- gdata::read.xls(xlsfile, sheet='Press Data', skip=1, header=T)[,1:5] # 5 cols
  mmd <- gdata::read.xls(xlsfile, sheet='Press Data (MinMax)', skip=1, header=T)[,1:5] # 5 cols
  temp <- gdata::read.xls(xlsfile, sheet='Temp Data', skip=1, header=T)[,1:5] # 5 cols
  temp_max <- gdata::read.xls(xlsfile, sheet='Temp Data (MinMax)', skip=1, header=T)[,1:5] # 5 cols
  
  ## MMD
  ## create max depth df for calc.bathy
  mmd$Date <- as.POSIXct(mmd$Date.Time, format='%m/%d/%y', tz='UTC')
  mmd$MaxDepth <- mmd$Max.Depth.m. * -1
  mmd <- mmd[,c('Date','MaxDepth')]
  mmd <- mmd[which(mmd$Date <= iniloc$date[2]),]
  udates <- unique(as.Date(mmd$Date))
  gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
  print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have MMD data...', sep=''))
  
  ## PDT
  ## create depth-temperature profile data for various 3d likelihood calcs
  pdt_dates <- seq.POSIXt(tag, pop, by = 'day') ## generate daily summary of depth-temp profiles
  depth$Date <- as.POSIXct(depth$Date.Time, format='%m/%d/%y %H:%M', tz='UTC')
  depth$Depth <- depth$Depth.m. * -1
  temp$Date <- as.POSIXct(temp$Date.Time, format='%m/%d/%y %H:%M', tz='UTC')
  pdt <- merge(depth, temp, by='Date')
  pdt <- pdt[,c('Date','Depth','Temp.C.')]
  names(pdt)[3] <- 'Temperature'
  ## generate depth-temp summary from time series
  pdt <- bin_TempTS(pdt, out_dates = pdt_dates, bin_res = 25)
  pdt <- pdt[,c('Date','Depth','MinTemp','MaxTemp')]
  udates <- unique(as.Date(pdt$Date))
  gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
  print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have PDT data...', sep=''))
  
  ## SST
  ## assume daily max temp is equivalent to SST
  sst <- temp_max
  sst$Date <- as.POSIXct(sst$Date.Time, format='%m/%d/%y', tz='UTC')
  sst$Temperature <- sst$Max.Temp.C.
  sst <- sst[,c('Date','Temperature')]
  udates <- unique(as.Date(sst$Date))
  gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
  print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have SST data...', sep=''))
  
  ## LIGHT LOCS
  #MWTxy = .getMWTxy(xlsfile, x0, xT, day0, dayT)
  MWTxy = locs[!is.na(locs[, 1]), 1:3]
  MWTxy$Date <- as.POSIXct(MWTxy$Date, format='%b %d, %Y', tz='UTC')
  #dateVec = seq(day0, dayT, by = 'day')
  didx = match(MWTxy$Date, dateVec)
  len = length(dateVec)
  MWTdata = as.data.frame(array(NA, c(len, 3)))
  didx = didx[!is.na(didx)]
  MWTdata[, 1] = dateVec
  MWTdata[, 2:3] = MWTxy[didx, 2:3]
  MWTdata[, 3] = -1 * (MWTdata[, 3])
  MWTdata[1, 2:3] = x0
  MWTdata[len, 2:3] = xT
  names(MWTdata) = c('Date','Latitude','Longitude')
  locs <- MWTdata
  locs$Error.Semi.minor.axis = .7 * 1000 * 111
  #locs$Error.Semi.minor.axis[didx] = fit$var.most.prob.track[,1]*1000*111
  #L.1 <- calc.gpe2(locs, locDates, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE, gpeOnly = TRUE)
  
  ## IT IS HIGHLY RECOMMENDED YOU CHECK AND MANUALLY FILTER OUTPUT LOCATIONS FROM MICROWAVE TAGS!!
  ## THESE ARE VERY RAW LIGHT-BASED ESTIMATES AND WHEN OUTLIERS ARE INCLUDED THEY CAN DRAMATICALLY IMPACT HMMOCE OUTPUTS
  ## here we started with 37 estimates and filtered to 8 reasonably useful ones
  plot(locs$Date, locs$Longitude, pch=16); world(add=T)
  if (aa == 1) locs <- locs[which(locs$Longitude < -40 & locs$Longitude > -80),]
  plot(locs$Longitude, locs$Latitude); world(add=T)
  udates <- unique(as.Date(locs$Date[which(!is.na(locs$Longitude))]))
  gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
  print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have Lightloc data...', sep=''))
  
  
  #invisible(readline(prompt=paste('Check plots. If passes QC, press [enter] to continue.', sep='')))
  
  
  if (aa == 1){
    sp.lim <- list(lonmin = -80,
                   lonmax = -60,
                   latmin = 15,
                   latmax = 45)
    
  } else if (aa == 2){
    sp.lim <- list(lonmin = -80,
                   lonmax = -60,
                   latmin = 15,
                   latmax = 45)
    
  } else if (aa == 3){
    sp.lim <- list(lonmin = -80,
                   lonmax = -60,
                   latmin = 15,
                   latmax = 45)
  }
    
  ## setup the spatial grid to base likelihoods on
  locs.grid <- setup.locs.grid(sp.lim, res='quarter')
  
  udates <- seq.Date(as.Date(tag), as.Date(pop), by = 'day')
  
  setwd('~/work/RCode/HMMoce_run/57508/')
  dir <- getwd()
  #load('./141259_tryHMMoce_20200724.rda')
  #load('./141259_L.res_20200727.rda')
  
  sst.dir <- paste0(dir, '/EnvData/sst/')
  if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
  for (i in 1:length(udates)) get.env(udates[i], filename='oisst', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)
  
  hycom.dir <- paste0(dir,'/EnvData/hycom/')
  #if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
  #for (i in 1:length(udates)) get.env(udates[i], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
  
  bathy.dir <- paste0(dir, '/EnvData/bathy/')
  if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
  #if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
  #bathy <- HMMoce::get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir, res=1)
  #} else{ ## OR (once downloaded and reading the .nc later)
  bathy <- irregular_ncToRaster(paste0(bathy.dir, 'bathy.nc'), varid = 'topo')
  #}
  
  
  L.light <- calc.lightloc(locs, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE)
  
  L.sst <- calc.sst(sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 1, focalDim = 3)
  
  
  # OCEAN HEAT CONTENT (INTEGRATED PDTS)
  L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  #L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
  
  # HYCOM PROFILE BASED LIKELIHOODS
  L.hycom.se <- calc.hycom(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
  #L.hycom.par <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:5], use.se = F, ncores=2)
  
  bathy_resamp <- raster::resample(bathy, L.sst) # or whatever grid makes sense to resample to
  L.bathy <- calc.bathy(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')
  
  ## make list of rasters
  L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.hycom.se = L.hycom.se, L.bathy = L.bathy)
  
  ## resample rasters
  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  save(L.res, file='57508_L.res_20200825.rda')
  save(L.rasters, file='57508_L.rasters_20200825.rda')
  
  load('57508_L.res_20200825.rda')
#}
  
  
  likVec <- c(1:5) 
  if (length(likVec) > 2){
    combine_idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F), utils::combn(likVec, 4, simplify=F))
  } else{
    combine_idx <- utils::combn(likVec, 2, simplify=F)
  }
  
  
  combine_idx <- combine_idx[!unlist(lapply(combine_idx, FUN=function(x) 3 %in% x & 4 %in% x))]
  
  # which of combine_idx combinations do you want to run?
  run.idx <- c(1:3,5,6,10:18)
  
  ## each run idx combination
  for (tt in run.idx[2:length(run.idx)]){
    
    L <- make.L(L.res$L.rasters[combine_idx[[tt]]], iniloc, dateVec)
    L.mle <- coarse.L(L, L.res$L.rasters)$L.mle
    g.mle <- coarse.L(L, L.res$L.rasters)$g.mle
    
    pars.optim <- opt.params(pars.init = c(2), 
                             lower.bounds = c(1), 
                             upper.bounds = c(8), 
                             g = L.res$g, 
                             #g = g.mle, 
                             L = L, 
                             #L = L.mle, 
                             alg.opt = 'optim', 
                             write.results = FALSE)
    
    #pars.optim <- opt.params(pars.init = c(2,2,1,1), 
    #                         lower.bounds = c(0.1, 0.1, 1, 1), 
    #                         upper.bounds = c(6, 6, 1, 1), 
    #                         g = L.res$g, 
    #                         #g = g.mle, 
    #                         L = L, 
    #                         #L = L.mle, 
    #                         alg.opt = 'optim', 
    #                         write.results = FALSE)
    ## about 22 mins per iteration on blue shark 141259
    ## about 1.5 mins with MLE grid
    ## final vals way different between the different grids
    
    #pars.nlminb <- opt.params(pars.init = c(2,.2,.6,.8), 
    #                          lower.bounds = c(0.1, 0.001, .1, .1), 
    #                          upper.bounds = c(5, .5, .9, .9), 
    #                          g = L.res$g, 
    #                          #g = g.mle, 
    #                          L = L, 
    #                          #L = L.mle, 
    #                          alg.opt = 'nlminb', 
    #                          write.results = FALSE)
    ## about 30 mins per iteration on blue shark 141259
    
    #ceiling(parallel::detectCores() * .9)
    #pars.ga <- opt.params(pars.init = c(2,2,1,1), 
    #                      lower.bounds = c(0.1, 0.1, 1, 1), 
    #                      upper.bounds = c(6, 6, 1, 1), 
    #                      g = L.res$g, 
    #                      #g = g.mle, 
    #                      L = L, 
    #                      #L = L.mle, 
    #                      alg.opt = 'ga', 
    #                      write.results = FALSE,
    #                      ncores = 2)
                          #ncores = ceiling(parallel::detectCores() * .9))
    ## about 7 hrs per iteration on blue shark 141259 w 4 cores
    ## about 1.6 hrs per iteration on blue shark 141259 w 15 cores
    ## about 2 mins with MLE grid but results way different
    #> pars.ga
    #$par
    #sigma1.ncell sigma2.ncell pswitch11 pswitch22
    #[1,]     4.937692    0.2290643 0.4023717  0.556778
    
    ## each pars approach
    #for (bb in c(1,3)){
    for (bb in c(1)){
      
      if (bb == 1){
        pars <- pars.optim$par
        #} else if (bb == 2){
        #   pars <- pars.nlminb$par
      } else if (bb == 3){
        pars <- pars.ga$par
      }
      
      for (ii in 1:2){
        
        #idx <- which(outVec_df[,1] == tt & outVec_df[,2] == bb & outVec_df[,3] == ii)
        #pars <- as.numeric(outVec_df[idx,4:7])
        
        if (length(pars) == 4){
          sigmas = pars[1:2]
          sizes = rep(ceiling(sigmas[1]*4),2)
          pb = pars[3:4]
          muadvs = c(0,0)
        } else if (length(pars) == 1){
          sigmas = pars[1]
          sizes = rep(ceiling(sigmas[1]*4),2)
          pb = NULL
          muadvs = c(0)
        }
        
        gausskern.PG <- ifelse(ii == 1, TRUE, FALSE)
        
        if(gausskern.PG){
          # behav 1
          if(sizes[1]%%2==0){sizes[1]=sizes[1]+1}
          ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
          while(ss<.999){
            sizes[1]=sizes[1]+2
            ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
          }
          K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
          rm(ss)
          K1=mask.K(K1)
          
          # behav 2
          if (!is.null(pb)){
            if(sizes[2]%%2==0){sizes[2]=sizes[2]+1}
            ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
            while(ss<.999){
              sizes[2]=sizes[2]+2
              ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
            }
            K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
            rm(ss)
            K2=mask.K(K2)
          }
          
          
        }else{
          sizes=rep(ceiling(sigmas[1]*4),2)
          K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
          if (!is.null(pb)) K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
        }
        
        
        if (!is.null(pb)){
          P <- matrix(c(pb[1],1-pb[1],1-pb[2],pb[2]),2,2,byrow=TRUE)
        } else{
          P <- NULL
        }
        
        # RUN THE FILTER STEP
        #f <- hmm.filter(L.res$g, L, K1, K2, maskL=FALSE, P)
        f <- hmm.filter(g=L.res$g, L=L, K=list(K1), maskL=FALSE, P=P, m=1)
        nllf <- -sum(log(f$psi[f$psi>0])) # negative log-likelihood
        
        # RUN THE SMOOTHING STEP
        s <- hmm.smoother(f, K=list(K1), L=L, P=P)
        
        # GET THE MOST PROBABLE TRACK
        tr <- calc.track(s, L.res$g, dateVec, iniloc)
        #plotHMM(s, tr, dateVec, ptt=paste(tt, bb, ii, sep='_'), save.plot = F)
        
        # WRITE OUT RESULTS
        if (is.null(pb)) pars <- c(pars, rep(NA, 3))
        outVec <- matrix(c(tt, bb, ii, pars, NLL = nllf), ncol = 8)
        res <- list(outVec = outVec, s = s, g = L.res$g, tr = tr, dateVec = dateVec, iniloc = iniloc)
        save(res, file=paste(paste(tt, bb, ii, sep = '_'), '-HMMoce_res_20200824_b1.rda', sep = ''))
        
        write.table(outVec, file=paste('141259_HMMoce_results_outVec_20200824.csv', sep = ''), sep = ',', append = T, row.names = F, col.names = F)
        
      } ## gausskern loop
      
    } ## par loop
    
  } ## run idx loop
  
  
  
  
  # set lik colors
  lik.breaks <- seq(0, 1, length.out = 25)
  lik.mid = lik.breaks[1:(length(lik.breaks)-1)]
  #lik.col = jet.colors(length(lik.breaks)-1) #[as.vector((dataT))]
  lik.col = terrain.colors(length(lik.breaks)-1, rev=T) #[as.vector((dataT))]
  
  f.sum <- apply(f$pred, 2:4, sum)
  s.sum <- apply(s, 2:4, sum)
  
  
  
  pdf(paste0(dir, '/', paste0(paste(tt, bb, ii, sep='_'), '_diagnose.pdf')),
      width=12, height=8, onefile=TRUE)
  
  #for (i in 1:10){
  for (i in 1:length(dateVec)){
    
    if (length(combine_idx[[tt]]) == 2) nf <- layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), nrow=6, ncol=2, byrow=F), widths=c(5,5), heights=c(4,4,4,4,4,4))
    if (length(combine_idx[[tt]]) == 3) nf <- layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow=F), widths=c(5,5), heights=c(4,4,4))
    if (length(combine_idx[[tt]]) == 4) nf <- layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,6,6,7,7), nrow=6, ncol=3, byrow=F), widths=c(5,5,5), heights=c(4,4,4,4,4,4))
    
    #layout.show(nf)
    
    #========
    ## plot likelihood 1
    #========
    par (mar=c(2,4,4,2))
    
    image(L.res$L.rasters[combine_idx[[tt]]][[1]][[i]], col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[combine_idx[[tt]][[1]]])
    axis(2); box()
    world(add=T)
    #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    
    
    #========
    ## plot likelihood 2
    #========
    par (mar=c(2,4,4,2))
    
    image(L.res$L.rasters[combine_idx[[tt]]][[2]][[i]], col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes = F, main=names(L.res$L.rasters)[combine_idx[[tt]][[2]]])
    axis(1); axis(2); box()
    world(add=T)
    #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    # legend
    #image(1, lik.mid, t(as.matrix(lik.mid)), breaks=lik.breaks, col=lik.col, axes=FALSE, xlab="",
    #      ylab=parse(text=paste('Temperature', "*degree~C", sep="")))
    #axis(2);box();
    
    if (length(combine_idx[[tt]]) >= 3){
      #========
      ## plot likelihood 3
      #========
      
      image(L.res$L.rasters[combine_idx[[tt]]][[3]][[i]], col = lik.col, breaks = lik.breaks,
            xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[combine_idx[[tt]][[3]]])
      axis(1); axis(2); box()
      world(add=T)
      #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
      
      
    }
    
    
    
    if (length(combine_idx[[tt]]) == 4){
      #========
      ## plot likelihood 4
      #========
      
      image(L.res$L.rasters[combine_idx[[tt]]][[4]][[i]], col = lik.col, breaks = lik.breaks,
            xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[combine_idx[[tt]][[4]]])
      axis(1); axis(2); box()
      world(add=T)
      #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
      
      
    }
    
    
    #========
    ## plot L
    #========
    image(L.res$g$lon[1,], L.res$g$lat[,1], L[i,,], col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes=F, main=dateVec[i])
    axis(1); axis(2); box()
    world(add=T)
    #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    #========
    ## plot filter
    #========
    image(L.res$g$lon[1,], L.res$g$lat[,1], f.sum[i,,] / max(f.sum[i,,]), col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes=F, main='filter')
    axis(1); axis(2); box()
    world(add=T)
    #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    
    
    #========
    ## plot smoother
    #========
    image(L.res$g$lon[1,], L.res$g$lat[,1], s.sum[i,,] / max(s.sum[i,,]), col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes=F, main='smoother')
    axis(1); axis(2); box()
    world(add=T)
    #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
  }
  dev.off()
  