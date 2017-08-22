# get the best track for each indiv
u.ptt <- unique(dat$ptt)
dat$final <- NA

for (i in 1:length(u.ptt)){
  ptt <- u.ptt[i]
  setwd(paste('~/ebs/Data/BaskingSharks/batch/', ptt, '/', sep=''))
  rev.idx <- which(dat$ptt == ptt)
  for (b in length(rev.idx):1){
    name.b <- dat$runname[rev.idx[b]]
    load(paste(name.b,'-HMMoce_res.rda', sep=''))
    
    # check iniloc
    if (any(!c(res$iniloc$lon[1] - 1 < round(res$tr$lon[1],0) & res$iniloc$lon[1] + 1 > round(res$tr$lon[1],0),
               res$iniloc$lat[1] - 1 < round(res$tr$lat[1],0) & res$iniloc$lat[1] + 1 > round(res$tr$lat[1],0),
               res$iniloc$lon[2] - 1 < round(res$tr$lon[nrow(res$tr)],0) & res$iniloc$lon[2] + 1 > round(res$tr$lon[nrow(res$tr)],0),
               res$iniloc$lat[2] - 1 < round(res$tr$lat[nrow(res$tr)],0) & res$iniloc$lat[2] + 1 > round(res$tr$lat[nrow(res$tr)],0)))){
      is.iniloc <- FALSE
    } else{
      is.iniloc <- TRUE
    }
    
    # check track distances
    # need to be less than 400 km between daily positions (based on 4 m/s sustained swimming 24 hrs)
    dists <- fields::rdist.earth.vec(res$tr[c(1:(nrow(res$tr)-1)),2:3], res$tr[c(2:nrow(res$tr)),2:3], miles=F)
    if (any(dists > 400)){
      is.trdist <- FALSE
    } else{
      is.trdist <- TRUE
    }
    
    if (is.iniloc & is.trdist){
      dat$final[rev.idx[b]] <- name.b
      break
    }
    
    if (b == 1) warning(paste('Nothing found for ', ptt, '...', sep=''))
    
  }
  
}

aws.s3::s3save(dat, bucket=bucketDir, object='all_bask_outvec_clean')

