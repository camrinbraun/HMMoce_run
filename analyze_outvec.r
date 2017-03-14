data <- read.table('~/Downloads/outVec_results.csv', sep=',', header=F, blank.lines.skip = F)
colnames(data) <- list('rownum','ptt', 'minBnd','migr.spd','rmselon.ti','rmselon.tib','rmselon.kf','rmselon.kfb','rmselon.gpe','rmselon.hmm',
                   'rmselat.ti','rmselat.tib','rmselat.kf','rmselat.kfb','rmselat.gpe','rmselat.hmm',
                   'gcdm.ti','gcdm.tib','gcdm.kf','gcdm.kfb','gcdm.gpe','gcdm.hmm',
                  'gcdsd.ti','gcdsd.tib','gcdsd.kf','gcdsd.kfb','gcdsd.gpe','gcdsd.hmm', 'L.idx',
                  'P.migr','P.resid')

#data <- data[which(data[,1] == 2),]
data <- data[,-grep('kf', names(data))]
data <- data[,-grep('ti', names(data))]
data <- data[which(data$minBnd == 10),]
data <- data[which(data$L.idx < 150),]
#elminate L.idx > 20 & < 50
data <- data[which(data$L.idx < 20 | data$L.idx > 50),]

library(dplyr)
summarise(group_by(data, minBnd))#, mgc=mean(gcdm.hmm), sdgc=sd(gcdm.hmm),
          mrm.lon=mean(rmselon.hmm), mrm.lat=mean(rmselat.hmm))
group_by(data, minBnd)

ggplot(data, aes(x=idx, y=gcdm.hmm, group=minBnd)) + geom_line(aes(col=as.factor(minBnd)))
ggplot(data, aes(gcdm.hmm, colour=as.factor(rownum))) + geom_density() +
  guides(colour=guide_legend(title='RUN (2=rescale)'))# + facet_wrap(~minBnd)
