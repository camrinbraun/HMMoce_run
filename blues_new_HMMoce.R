devtools::load_all('~/ebs/RCode/HMMoce')
library(raster)
library(foieGras)
library(tidyverse)

fList <- list.files('~/ebs/Data/data_org/', full.names = T, recursive = T)
fList <- fList[grep('141259', fList)]

meta <- read.table('../nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, stringsAsFactors = F)
meta[which(meta$uid_no == 9),]

rds_files <- list.files('~/ebs/RCode/nip_drake/move_res/', full.names = TRUE)
rds_ids <- rep(NA, length.out=length(rds_files))
for (i in 1:length(rds_ids)) rds_ids[i] <- substr(rds_files[i], stringr::str_locate(rds_files[i], '//')[2] + 1, stringr::str_locate(rds_files[i], '_move')[1] - 1)
idx <- which(rds_ids == '160424_2015_141261')

devtools::load_all('~/ebs/RCode/analyzePSAT')
move_eddies <- readRDS(rds_files[idx])
pred <- move_eddies$pred
raw <- move_eddies$raw

iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27,
                              10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
names(iniloc) <- list('day','month','year','lat','lon')
iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
tag <- iniloc$date[1]
pop <- iniloc$date[2]

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- seq.POSIXt(tag, pop, by = '24 hours')


sstFile <- system.file("extdata", "141259-SST.csv", package = "HMMoce")
#sstFile <- fList[grep('SST', fList)]
tag.sst <- read.wc(sstFile, type = 'sst', tag=tag, pop=pop, verbose=F)
tag.sst <- tag.sst[,c('Date','Depth','Temperature')]
head(tag.sst)


pdtFile <- system.file("extdata", "141259-PDTs.csv", package = "HMMoce")
#pdtFile <- fList[grep('PDTs', fList)]
pdt <- read.wc(pdtFile, type = 'pdt', tag = tag, pop = pop, verbose = F)
pdt <- pdt[,c('Date','Depth','MinTemp','MaxTemp')]
head(pdt)

llFile <- system.file("extdata", "141259-Locations-GPE2.csv", package = "HMMoce")
lightloc <- read.table(llFile, sep = ',', header = T, blank.lines.skip = F)
lightloc <- lightloc[which(lightloc$Type != 'Argos'),]
lightloc <- lightloc[,c('Date','Longitude','Error.Semi.minor.axis','Latitude','Error.Semi.major.axis','Offset','Offset.orientation')]
lightloc$Date <- as.POSIXct(lightloc$Date, format = findDateFormat(lightloc$Date))
head(lightloc)

mmdFile <- system.file("extdata", "141259-MinMaxDepth.csv", package = "HMMoce")
mmd <- read.table(mmdFile, sep = ',', header = T, blank.lines.skip = F)[,c('Date','MinDepth','MaxDepth')]
mmd$Date <- as.POSIXct(mmd$Date, format = findDateFormat(mmd$Date))
head(mmd)


sp.lim <- list(lonmin = -75,
               lonmax = -30,
               latmin = 20,
               latmax = 47)

## setup the spatial grid to base likelihoods on
locs.grid <- setup.locs.grid(sp.lim, res='quarter')

## we only need udates here if dateVec resolution is finer than 1 day as we typically only download env data as fine as daily resolution
udates <- seq.Date(as.Date(tag), as.Date(pop), by = 'day')

setwd('~/ebs/Data/HMMoce_run/141259/')
dir <- getwd()
#load('./141259_tryHMMoce_20200724.rda')
load('./141259_L.res_20200727.rda')

sst.dir <- paste0(dir, '/EnvData/sst/')
if (!dir.exists(sst.dir)) dir.create(sst.dir, recursive = TRUE)
for (i in 1:length(udates)) get.env(udates[i], filename='oisst', type = 'sst', sst.type='oi', spatLim = sp.lim, save.dir = sst.dir)

hycom.dir <- paste0(dir,'/EnvData/hycom/')
if (!dir.exists(hycom.dir)) dir.create(hycom.dir, recursive = TRUE)
for (i in 1:length(udates)) get.env(udates[i], filename='hycom', type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

bathy.dir <- paste0(dir, '/EnvData/bathy/')
if (!dir.exists(bathy.dir)) dir.create(bathy.dir, recursive = TRUE)
#if (!file.exists(paste0(bathy.dir, 'bathy.nc'))){
bathy <- HMMoce::get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, folder = bathy.dir, res=1)
#} else{ ## OR (once downloaded and reading the .nc later)
  bathy <- irregular_ncToRaster(paste0(bathy.dir, 'bathy.nc'), varid = 'topo')
#}


L.light <- calc.lightloc(lightloc, locs.grid = locs.grid, dateVec = dateVec, errEll = FALSE)

L.sst <- calc.sst.par(tag.sst, filename='oisst', sst.dir = sst.dir, dateVec = dateVec, sens.err = 3, focalDim = 3)


# OCEAN HEAT CONTENT (INTEGRATED PDTS)
L.ohc <- calc.ohc(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)
#L.ohc <- calc.ohc.par(pdt, filename='hycom', ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = F)

# WORLD OCEAN ATLAS-BASED LIKELIHOODS
#L.woa <- calc.woa(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)
#L.woa <- calc.woa.par(pdt, woa.data = woa.quarter, sp.lim=sp.lim, focalDim = 9, dateVec = dateVec, use.se = T)

# HYCOM PROFILE BASED LIKELIHOODS
L.hycom.se <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec, use.se = T)
#L.hycom.par <- calc.hycom.par(pdt, filename='hycom', hycom.dir, focalDim = 9, dateVec = dateVec[1:5], use.se = F, ncores=2)

#L.bt <- calc.bottomTemp(tag.bt, dateVec[1:5], focalDim = 3, sens.err = 1, bt.dir = sst.dir, filename = 'oisst', varName = 'sst')

bathy_resamp <- raster::resample(bathy, L.sst) # or whatever grid makes sense to resample to
L.bathy <- calc.bathy(mmd, bathy_resamp, dateVec, focalDim = 3, sens.err = 5, lik.type = 'max')

## make list of rasters
L.rasters <- list(L.light = L.light, L.sst = L.sst, L.ohc = L.ohc, L.hycom.se = L.hycom.se, L.bathy = L.bathy)

## resample rasters
resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
save(L.res, file='141259_L.res_20200727.rda')
save(L.rasters, file='141259_L.rasters_20200727.rda')


likVec <- c(1:5) 
if (length(likVec) > 2){
  combine_idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F), utils::combn(likVec, 4, simplify=F))
} else{
  combine_idx <- utils::combn(likVec, 2, simplify=F)
}


combine_idx <- combine_idx[!unlist(lapply(combine_idx, FUN=function(x) 3 %in% x & 4 %in% x))]

# which of combine_idx combinations do you want to run?
run.idx <- c(1:3,5,6,10:18)

outVec_df <- read.table(paste('141259_HMMoce_results_outVec_20200824.csv', sep=''), sep=',', header=F)
## each run idx combination
for (tt in run.idx[6:length(run.idx)]){
#for (tt in run.idx[5:4]){
    
  L <- make.L(L.res$L.rasters[combine_idx[[tt]]], iniloc, dateVec)
  L.mle <- coarse.L(L, L.res$L.rasters)$L.mle
  g.mle <- coarse.L(L, L.res$L.rasters)$g.mle
}

  pars.optim <- opt.params(pars.init = c(2,.2,.6,.8), 
                           lower.bounds = c(0.1, 0.001, .1, .1),  
                           upper.bounds = c(6, .6, .9, .9), 
                           g = L.res$g, 
                           #g = g.mle, 
                           L = L, 
                           #L = L.mle, 
                           alg.opt = 'optim', 
                           write.results = FALSE)
  
  pars.optim.mle <- opt.params(pars.init = c(2,.2,.6,.8), 
                           lower.bounds = c(0.1, 0.001, .1, .1),  
                           upper.bounds = c(6, .6, .9, .9), 
                           #g = L.res$g, 
                           g = g.mle, 
                           #L = L, 
                           L = L.mle, 
                           alg.opt = 'optim', 
                           write.results = FALSE)
  
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
  pars.ga.one <- opt.params(pars.init = c(2), 
                        lower.bounds = c(1), 
                        upper.bounds = c(8), 
                        g = L.res$g, 
                        #g = g.mle, 
                        L = L, 
                        #L = L.mle, 
                        alg.opt = 'ga', 
                        write.results = FALSE,
                        ncores = ceiling(parallel::detectCores() * .9))
                        #ncores = 2)
  
  pars.ga.mle <- opt.params(pars.init = c(2,.2,.6,.8), 
                        lower.bounds = c(0.1, 0.001, .1, .1), 
                        upper.bounds = c(6, .6, .9, .9), 
                        #g = L.res$g, 
                        g = g.mle, 
                        #L = L, 
                        L = L.mle, 
                        alg.opt = 'ga', 
                        write.results = FALSE,
                        ncores = ceiling(parallel::detectCores() * .9))
                        #ncores = 2)

  #pars.ga <- opt.params(pars.init = c(2,.2,.6,.8), 
  #                          lower.bounds = c(0.1, 0.001, .1, .1), 
  #                          upper.bounds = c(6, .6, .9, .9), 
  #                          g = L.res$g, 
                            #g = g.mle, 
  #                          L = L, 
                            #L = L.mle, 
  #                          alg.opt = 'ga', 
  #                          write.results = FALSE,
  #                          ncores = ceiling(parallel::detectCores() * .9))

  ## about 7 hrs per iteration on blue shark 141259 w 4 cores
  ## about 1.6 hrs per iteration on blue shark 141259 w 15 cores
  ## about 2 mins with MLE grid but results way different
  #> pars.ga
  #$par
  #sigma1.ncell sigma2.ncell pswitch11 pswitch22
  #[1,]     4.937692    0.2290643 0.4023717  0.556778
  
  
  #pars.optim; pars.nlminb; pars.ga
for (tt in run.idx[5:length(run.idx)]){
    #for (tt in run.idx[5:4]){
    
    L <- make.L(L.res$L.rasters[combine_idx[[tt]]], iniloc, dateVec)
    #L.mle <- coarse.L(L, L.res$L.rasters)$L.mle
    #g.mle <- coarse.L(L, L.res$L.rasters)$g.mle
    
  ## each pars approach
  for (bb in c(3)){
  #for (bb in c(3)){
      
    if (bb == 1){
      pars <- pars.optim$par
    } else if (bb == 2){
      pars <- pars.optim.mle$par
    } else if (bb == 3){
      
      pars.ga <- opt.params(pars.init = c(2,.2,.6,.8), 
                                lower.bounds = c(0.1, 0.001, .1, .1), 
                                upper.bounds = c(6, .6, .9, .9), 
                                g = L.res$g, 
      #g = g.mle, 
                                L = L, 
      #L = L.mle, 
                                alg.opt = 'ga', 
                                write.results = FALSE,
                                ncores = ceiling(parallel::detectCores() * .9))
      
      pars <- pars.ga$par
    } else if (bb == 4){
      pars <- pars.ga.mle$par
    } else if (bb == 5){
      pars <- pars.ga.one$par
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
      
      ## set transition matrix, if applicable
      if (!is.null(pb)){
        P <- matrix(c(pb[1],1-pb[1],1-pb[2],pb[2]),2,2,byrow=TRUE)
      } else{
        P <- NULL
      }
      
      # RUN THE FILTER STEP
      if (!is.null(pb)){
        K=list(K1,K2)
        f <- hmm.filter(g=L.res$g, L=L, K=K, maskL=FALSE, P=P, m=2)
      } else{
        K=list(K1)
        f <- hmm.filter(g=L.res$g, L=L, K=K, maskL=FALSE, P=P, m=1)
      }
      nllf <- -sum(log(f$psi[f$psi>0])) # negative log-likelihood
      
      # RUN THE SMOOTHING STEP
      s <- hmm.smoother(f, K=K, L=L, P=P)
      
      # GET THE MOST PROBABLE TRACK
      tr <- try(calc.track(s, g=L.res$g, dateVec, iniloc), TRUE)
      if (class(tr) == 'try-error') next
      #plotHMM(s, tr, dateVec, ptt=paste(tt, bb, ii, sep='_'), save.plot = F)
      
      # WRITE OUT RESULTS
      if (is.null(pb)) pars <- c(pars, rep(NA, 3))
      outVec <- matrix(c(tt, bb, ii, pars, NLL = nllf), ncol=8)
      res <- list(outVec = outVec, s = s, g = L.res$g, tr = tr, dateVec = dateVec, iniloc = iniloc)
      save(res, file=paste(paste(tt, bb, ii, sep='_'), '-HMMoce_res_20200828.rda', sep=''))

      write.table(outVec, file=paste('141259_HMMoce_results_outVec_20200824.csv', sep=''), sep=',', append=T, row.names = F, col.names = F)
      
    } ## gausskern loop
    
  } ## par loop
  
} ## run idx loop



raw <- raw %>% filter(date <= iniloc$date[2])

raw[nrow(raw) + 1,] <- c('160424_2015_141261.1', NA, '3', iniloc$lon[2], iniloc$lat[2])
raw$date[nrow(raw)] <- iniloc$date[2]
raw$lon <- as.numeric(raw$lon); raw$lat <- as.numeric(raw$lat)

df.locs <- split(raw, raw$id)
ssm_fit <- lapply(df.locs, function(x) foieGras::fit_ssm(x, model='crw', time.step = 24, vmax=10))

## grab predicted locations output from fit_ssm
plocs <- lapply(ssm_fit, FUN=function(x) foieGras::grab(x, what = "p", as_sf = FALSE)) %>%
  do.call(rbind, .) %>%
  tbl_df() %>%
  mutate(id = as.character(id)) %>% group_by(id)

## subsample from predicted locations, if applicable
res_out <- 24
if (!is.null(res_out)){
  plocs <- split(plocs, plocs$id)
  
  for (i in 1:length(plocs)){
    res.i <- as.numeric(difftime(plocs[[i]]$date[2], plocs[[i]]$date[1], units = 'hours'))
    res_slice <- res_out / res.i
    if (res_slice < 1) res_slice <- 1
    plocs[[i]] <- plocs[[i]] %>% do(slice(., seq(1, n(), by = res_slice)))
  }
  
  plocs <- plocs %>%
    do.call(rbind, .) %>%
    tbl_df() %>%
    mutate(id = as.character(id)) %>% group_by(id)
  
}

mpm_fit <- try(fit_mpm(plocs[,c('id','date','lon','lat')], model = 'mpm'), silent=TRUE)
plocs$g <- NA
for (tt in 1:length(mpm_fit$mpm)){
  
  if (mpm_fit$converged[tt]){
    temp_g <- mpm_fit$mpm[[tt]]$fitted$g
    plocs$g[which(plocs$id == unique(plocs$id)[tt])] <- temp_g
    
  } else {
    #temp_g <- rep(NA, length.out = nrow(plocs[which(plocs$id == unique(plocs$id)[tt]),]))
  }
  
}

## convert to same scale as hmmoce behav estimate "p"
plocs$p_equiv <- plocs$g * -1 + 1
plocs$dateVec <- findInterval(plocs$date, dateVec)


out <- read.table('141259_HMMoce_results_outVec_20200824.csv', sep=',', header=F)
names(out) <- c('likvec','pars','ii','pars1','pars2','pars3','pars4','nll')
out[c((ncol(out) + 1):(ncol(out) + 8))] <- NA
names(out)[9:ncol(out)] <- c('rmse.lon','rmse.lat','gcd_mean','gcd_sd','gcd_median','gcd_min','gcd_max','rmse.behav')
out$likvec <- as.character(out$likvec)

for (i in 1:nrow(out)){
  #stringr::str_locate(out$likvec[i], '_')
  if (nchar(out$likvec[i]) > 2){
    result <- try(load(paste0(substr(out$likvec[i],1,1), substr(out$likvec[i],3,5), '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200828.rda')), silent=TRUE)
    if (class(result) == 'try-error'){
      result <- try(load(paste0(substr(out$likvec[i],1,1), substr(out$likvec[i],3,5), '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200824.rda')), silent=TRUE)
      if (class(result) == 'try-error'){
        result <- try(load(paste0(substr(out$likvec[i],1,1), substr(out$likvec[i],3,5), '_', out$pars[i], '_', out$ii[i], '-HMMoce_res.rda')), silent=TRUE)
      }
    }
  } else{
    result <- try(load(paste0(out$likvec[i], '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200828.rda')), silent=TRUE)
    if (class(result) == 'try-error'){
      result <- try(load(paste0(out$likvec[i], '_', out$pars[i], '_', out$ii[i], '-HMMoce_res_20200824.rda')), silent=TRUE)
      if (class(result) == 'try-error'){
        result <- try(load(paste0(out$likvec[i], '_', out$pars[i], '_', out$ii[i], '-HMMoce_res.rda')), silent=TRUE)
      }
    }
  } 
  
  if (class(result) == 'try-error') next
  # res <- list(outVec = outVec, s = s, g = L.res$g, tr = tr, dateVec = dateVec, iniloc = iniloc)

  comp <- compareTracks(res$tr, plocs, dateVec, plot=FALSE)
  
  out[i,9:ncol(out)] <- unlist(comp)
  print(i)
  rm(res)
}

summary(out)
write.table(out, '141259_HMMoce_results_outVec_20200824_compare.csv', sep=',', col.names=T, row.names=F)

out %>%
  select(-gcd_min, -gcd_max) %>%
  gather(-likvec, -pars, -ii, -pars1, -pars2, -pars3, -pars4, -nll, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = nll, color = factor(likvec), shape = factor(pars), size=ii*.5)) +
  geom_point() + ylim(90,200) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

out %>%
  select(pars, rmse.lon) %>%
  filter(pars %in% c(1,3,5)) %>%
  #gather(-likvec, -pars, -ii, -pars1, -pars2, -pars3, -pars4, -nll, key = "var", value = "value") %>% 
  ggplot(aes(as.factor(pars), rmse.lon)) +
  geom_boxplot() 

out_ii <- out %>%
  group_by(likvec, pars) %>% summarise(nll1 = nll[which(ii==1)], nll2 = nll[which(ii==2)],
                                       gcd_mean1 = gcd_mean[which(ii==1)], gcd_mean2 = gcd_mean[which(ii==2)])
out_ii$diff <- out_ii$gcd_mean2 - out_ii$gcd_mean1  
out_ii$diff_nll <- out_ii$nll2 - out_ii$nll1  

select(ii, rmse.lon) %>%
  #filter(pars %in% c(1,3,5)) %>%
  #gather(-likvec, -pars, -ii, -pars1, -pars2, -pars3, -pars4, -nll, key = "var", value = "value") %>% 
  ggplot(aes(as.factor(ii), rmse.lon)) +
  geom_boxplot() 

out %>%
  select(likvec, rmse.lon) %>%
  #gather(-likvec, -pars, -ii, -pars1, -pars2, -pars3, -pars4, -nll, key = "var", value = "value") %>% 
  ggplot(aes(as.factor(likvec), rmse.lon)) +
  geom_boxplot()

for (i in 1:length(dateVec)){
  image.plot(res$g$lon[1,], res$g$lat[,1], res$s[1,i,,])
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  invisible(readline(prompt=paste('Check plots. If passes QC, press [enter] to continue.', sep='')))
  
}




setwd('~/ebs/Data/HMMoce_run/141259/')
dir <- getwd()
#load('./141259_tryHMMoce_20200724.rda')
#save(L.res, file='141259_L.res_20200724.rda')
#save(plocs, file='141259_plocs_20200724.rda')

load('141259_L.res_20200724.rda')
load('141259_plocs_20200724.rda')

# set lik colors
lik.breaks <- seq(0, 1, length.out = 25)
lik.mid = lik.breaks[1:(length(lik.breaks)-1)]
#lik.col = jet.colors(length(lik.breaks)-1) #[as.vector((dataT))]
lik.col = terrain.colors(length(lik.breaks)-1, rev=T) #[as.vector((dataT))]

#png(paste0(dir, '/141259_diagnose.png'),
#    width=5000, height=3000, res=300, onefile=TRUE)
tt=12; bb=1; ii=1
load(paste(paste(tt, bb, ii, sep='_'), '-HMMoce_res.rda', sep=''))
iniloc <- res$iniloc; dateVec <- res$dateVec
tr <- res$tr

likVec <- c(1:5) 
if (length(likVec) > 2){
  combine_idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
} else{
  combine_idx <- utils::combn(likVec, 2, simplify=F)
}


L <- make.L(L.res$L.rasters[combine_idx[[tt]]], iniloc, dateVec)

pars <- res$outVec[4:7]
sigmas=pars[1:2]
sizes=ceiling(sigmas*4)
pb=pars[3:4]
muadvs=c(0,0)

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
    if(sizes[2]%%2==0){sizes[2]=sizes[2]+1}
    ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
    while(ss<.999){
      sizes[2]=sizes[2]+2
      ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
    }
    K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
    rm(ss)
    K2=mask.K(K2)
    
  }else{
    sizes=rep(ceiling(sigmas[1]*4),2)
    K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
    K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
  }
  
  
  P <- matrix(c(pb[1],1-pb[1],1-pb[2],pb[2]),2,2,byrow=TRUE)
  

f <- hmm.filter(L.res$g, L, K1, K2, maskL=FALSE, P)

f.sum <- apply(f$pred, 2:4, sum)
s.sum <- apply(res$s, 2:4, sum)


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
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  

  
  #========
  ## plot likelihood 2
  #========
  par (mar=c(2,4,4,2))
  
  image(L.res$L.rasters[combine_idx[[tt]]][[2]][[i]], col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes = F, main=names(L.res$L.rasters)[combine_idx[[tt]][[2]]])
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
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
    points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    
  }
  
  
  
  if (length(combine_idx[[tt]]) == 4){
    #========
    ## plot likelihood 4
    #========
    
    image(L.res$L.rasters[combine_idx[[tt]]][[4]][[i]], col = lik.col, breaks = lik.breaks,
          xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[combine_idx[[tt]][[4]]])
    axis(1); axis(2); box()
    world(add=T)
    points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
    
    
  }
  
  
  #========
  ## plot L
  #========
  image(L.res$g$lon[1,], L.res$g$lat[,1], L[i,,], col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes=F, main=dateVec[i])
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
  #========
  ## plot filter
  #========
  image(L.res$g$lon[1,], L.res$g$lat[,1], f.sum[i,,] / max(f.sum[i,,]), col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes=F, main='filter')
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
  
  
  #========
  ## plot smoother
  #========
  image(L.res$g$lon[1,], L.res$g$lat[,1], s.sum[i,,] / max(s.sum[i,,]), col = lik.col, breaks = lik.breaks,
        xlab='', ylab='', axes=F, main='smoother')
  axis(1); axis(2); box()
  world(add=T)
  points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
  
}
dev.off()
