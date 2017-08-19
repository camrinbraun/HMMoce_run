# writes out some stats on data coverage
read.wc2 <- function(ptt, wd = getwd(), tag, pop, type = 'sst'){
  
  if(substr(wd, nchar(wd), nchar(wd)) == '/'){
    
  } else{
    wd <- paste(wd, '/', sep='')
  }
  
  if(type == 'pdt'){
    # READ IN PDT DATA FROM WC FILES
    data <- utils::read.table(paste(wd, ptt,'-PDTs.csv', sep=''), sep=',',header=T,blank.lines.skip=F, skip = 0)
    print(paste('If read.wc() fails for type=pdt, check the number of column headers in the PDTs.csv file.'))
    data <- HMMoce:::extract.pdt(data)
    dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- dts >= (tag + d1) & dts <= (pop - d1)
    data <- data[didx,]; dts <- dts[didx]
    
    # check for days with not enough data for locfit
    data1 <- data
    data1$dts <- as.Date(dts)
    dt.cut <- data.frame(group_by(data1, dts) %>% summarise(n=n()))
    dt.cut <- dt.cut[which(dt.cut[,2] < 3),1]
    if (length(dt.cut) == 0){
      
    } else{
      data <- data1[-which(data1$dts %in% dt.cut),]
    }
    
    # get unique dates
    udates <- unique(as.Date(data$Date))
    
    # get data gaps
    print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have PDT data...', sep=''))
    data.count <- length(which(as.Date(seq(tag, pop, 'day')) %in% udates))
    deploydur <- length(seq(tag, pop, 'day'))
    gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
    print(paste('Data gaps are ', paste(gaps[gaps > 1], collapse=', '), ' days in PDT...'))
    gaps <- gaps[gaps > 1]
    
  } else if(type == 'sst'){
    # READ IN TAG SST FROM WC FILES
    data <- utils::read.table(paste(wd, ptt, '-SST.csv', sep=''), sep=',',header=T, blank.lines.skip=F)
    dts <- as.POSIXct(data$Date, format = findDateFormat(data$Date))
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- dts >= (tag + d1) & dts <= (pop - d1)
    data <- data[didx,]
    if (length(data[,1]) <= 1){
      stop('Something wrong with reading and formatting of tags SST data. Check date format.')
    }
    dts <- as.Date(as.POSIXct(data$Date, format = findDateFormat(data$Date)))
    udates <- unique(dts)
    data <- data[,c(1:11)]
    
    # get data gaps
    print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have GPE data...', sep=''))
    data.count <- length(which(as.Date(seq(tag, pop, 'day')) %in% udates))
    deploydur <- length(seq(tag, pop, 'day'))
    gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
    print(paste('Data gaps are ', paste(gaps[gaps > 1], collapse=', '), ' days...'))
    gaps <- gaps[gaps > 1]
    
  } else if(type == 'light'){
    # READ IN LIGHT DATA FROM WC FILES
    data <- utils::read.table(paste(wd,'/', ptt, '-LightLoc.csv', sep=''), sep=',',header=T, blank.lines.skip=F,skip=2)
    if(!any(grep('depth', names(data), ignore.case=T))) data <- utils::read.table(paste(wd,'/', ptt, '-LightLoc.csv', sep=''), sep=',',header=T, blank.lines.skip=F,skip=1)
    if(!any(grep('depth', names(data), ignore.case=T))) data <- utils::read.table(paste(wd,'/', ptt, '-LightLoc.csv', sep=''), sep=',',header=T, blank.lines.skip=F,skip=0)
    data <- data[which(!is.na(data[,1])),]
    
    #dts <- as.POSIXct(data$Day, format = findDateFormat(data$Day), tz = 'UTC')
    dts <- as.POSIXct(data$Day, format = '%d-%b-%y', tz = 'UTC')
    
    if(as.Date(dts[1]) > as.Date(Sys.Date()) | as.Date(dts[1]) < '1990-01-01'){
      stop('Error: dates not parsed correctly.')    
    }
    
    d1 <- as.POSIXct('1900-01-02') - as.POSIXct('1900-01-01')
    didx <- (dts > (tag + d1)) & (dts < (pop - d1))
    data <- data[didx,]
    data$dts <- as.Date(dts[didx])
    udates <- unique(as.Date(dts))
    
    # get data gaps
    print(paste(length(which(as.Date(seq(tag, pop, 'day')) %in% udates)), ' of ', length(seq(tag, pop, 'day')), ' deployment days have GPE data...', sep=''))
    gaps <- diff(c(as.Date(tag), udates, as.Date(pop)), units='days')
    print(paste('Data gaps are ', paste(gaps[gaps > 1], collapse=', '), ' days...'))
    
  }
  
  return(list(data = data, udates = udates, data.count=data.count, deploydur=deploydur, gaps=gaps))
  
}

