for (i in 4:nrow(meta)){

  work_dir <- paste0(base_dir,'/', meta$instrument_name[i], '/')
  setwd(work_dir)
  
  ## tracks not ran
  if (length(grep('1190241', meta$instrument_name[i])) > 0) next
  if (length(grep('1090269', meta$instrument_name[i])) > 0) next
  
  fList <- list.files()
  
  ## skip if compare pdf exists
  if (length(unlist(stringr::str_locate_all(fList, '_compare.pdf'))) > 0) next
  
  load(file=fList[grep('L.res', fList)])

  ## temporal bounds
  iniloc <- data.frame(matrix(c(lubridate::day(meta$time_coverage_start[i]),
                                lubridate::month(meta$time_coverage_start[i]),
                                lubridate::year(meta$time_coverage_start[i]),
                                meta$geospatial_lat_start[i],
                                make360(meta$geospatial_lon_start[i]),
                                lubridate::day(meta$time_coverage_end[i]),
                                lubridate::month(meta$time_coverage_end[i]),
                                lubridate::year(meta$time_coverage_end[i]),
                                meta$geospatial_lat_end[i],
                                make360(meta$geospatial_lon_end[i])),
                              nrow = 2, ncol = 5, byrow = T))
  
  names(iniloc) <- list('day','month','year','lat','lon')
  iniloc$date <- as.POSIXct(paste(iniloc$year, iniloc$month, iniloc$day, sep='-'), tz='UTC')
  tag <- iniloc$date[1]
  pop <- iniloc$date[2]
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  
  ## cut the sst and ohc data from a subset of individual's likelihoods due to external thermistor failure
  if (length(grep('B2381', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2005-05-11', tz='UTC')
  if (length(grep('D1464', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2007-04-29', tz='UTC')
  if (length(grep('A1973', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2004-08-14', tz='UTC')
  if (length(grep('A2088', meta$instrument_name[i])) > 0) cut_date <- as.POSIXct('2004-12-01', tz='UTC')
  
  if (length(grep('390191', meta$instrument_name[i])) > 0){
    cut_date <- as.POSIXct('2004-01-19', tz='UTC')
    iniloc$day[2] <- 19
    iniloc$month[2] <- 01
    iniloc$lat[2] <- 24.248
    iniloc$lon[2] <- 243.64
    iniloc$date[2] <- cut_date
    pop <- cut_date
    dateVec <- seq.POSIXt(tag, pop, by = '24 hours')
  }
  
  if (exists('cut_date')) dateVec <- seq.POSIXt(tag, cut_date, by = '24 hours')
  
  
  run_idx <- list(c(1:4),
                  c(1,2,4),
                  c(1,3,4))
  
  for (tt in 1){
    #load(file=paste0(meta$instrument_name[i],'_L_tt', tt, '_20201221.rda'))
    load(file=fList[max(grep('L_tt1_2', fList))])
    
    # set lik colors
    lik.breaks <- seq(0, 1, length.out = 25)
    lik.mid = lik.breaks[1:(length(lik.breaks)-1)]
    #lik.col = jet.colors(length(lik.breaks)-1) #[as.vector((dataT))]
    lik.col = terrain.colors(length(lik.breaks)-1, rev=T) #[as.vector((dataT))]
    
    pdf(paste0(meta$instrument_name[i], '_L_tt', tt, '_compare2.pdf'),
        width=12, height=8, onefile=TRUE)
    
    #for (i in 1:10){
    for (b in 1:length(dateVec)){
      
      if (length(run_idx[[tt]]) == 2) nf <- layout(matrix(c(1,1,1,2,2,2,3,3,4,4,5,5), nrow=6, ncol=2, byrow=F), widths=c(5,5), heights=c(4,4,4,4,4,4))
      if (length(run_idx[[tt]]) == 3) nf <- layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow=F), widths=c(5,5), heights=c(4,4,4))
      if (length(run_idx[[tt]]) == 4) nf <- layout(matrix(c(1,2,3,4,5,5), nrow=3, ncol=2, byrow=F), widths=c(9,9), heights=c(4,4,4))
      
      #layout.show(nf)
      
      #========
      ## plot likelihood 1
      #========
      par (mar=c(2,4,4,2))
      
      raster::image(L.res$L.rasters[run_idx[[tt]]][[1]][[b]], col = lik.col, breaks = lik.breaks,
                    xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[run_idx[[tt]][[1]]])
      axis(2); box()
      world(add=T, wrap=c(0,360))
      #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
      
      
      
      #========
      ## plot likelihood 2
      #========
      par (mar=c(2,4,4,2))
      
      raster::image(L.res$L.rasters[run_idx[[tt]]][[2]][[b]], col = lik.col, breaks = lik.breaks,
                    xlab='', ylab='', axes = F, main=names(L.res$L.rasters)[run_idx[[tt]][[2]]])
      axis(1); axis(2); box()
      world(add=T, wrap=c(0,360))
      #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
      
      # legend
      #image(1, lik.mid, t(as.matrix(lik.mid)), breaks=lik.breaks, col=lik.col, axes=FALSE, xlab="",
      #      ylab=parse(text=paste('Temperature', "*degree~C", sep="")))
      #axis(2);box();
      
      if (length(run_idx[[tt]]) >= 3){
        #========
        ## plot likelihood 3
        #========
        
        raster::image(L.res$L.rasters[run_idx[[tt]]][[3]][[b]], col = lik.col, breaks = lik.breaks,
                      xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[run_idx[[tt]][[3]]])
        axis(1); axis(2); box()
        world(add=T, wrap=c(0,360))
        #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
        
        
      }
      
      
      
      if (length(run_idx[[tt]]) == 4){
        #========
        ## plot likelihood 4
        #========
        
        raster::image(L.res$L.rasters[run_idx[[tt]]][[4]][[b]], col = lik.col, breaks = lik.breaks,
                      xlab='', ylab='', axes=F, main=names(L.res$L.rasters)[run_idx[[tt]][[4]]])
        axis(1); axis(2); box()
        world(add=T, wrap=c(0,360))
        #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
        
        
      }
      
      
      #========
      ## plot L
      #========
      raster::image(L.res$g$lon[1,], L.res$g$lat[,1], L[b,,], col = lik.col, breaks = lik.breaks,
                    xlab='', ylab='', axes=F, main=dateVec[b])
      axis(1); axis(2); box()
      world(add=T, wrap=c(0,360))
      #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
      
      #========
      ## plot filter
      #========
      #image(L.res$g$lon[1,], L.res$g$lat[,1], f.sum[i,,] / max(f.sum[i,,]), col = lik.col, breaks = lik.breaks,
      #      xlab='', ylab='', axes=F, main='filter')
      #axis(1); axis(2); box()
      #world(add=T, wrap=c(0,360))
      #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
      
      
      
      #========
      ## plot smoother
      #========
      #image(L.res$g$lon[1,], L.res$g$lat[,1], s.sum[i,,] / max(s.sum[i,,]), col = lik.col, breaks = lik.breaks,
      #      xlab='', ylab='', axes=F, main='smoother')
      #axis(1); axis(2); box()
      #world(add=T, wrap=c(0,360))
      #points(plocs$lon[which(plocs$dateVec == i)], plocs$lat[which(plocs$dateVec == i)])
      
    }
    dev.off()
    
  }
  
  rm(cut_date); rm(L.res); gc()
}
