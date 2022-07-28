data <- read.csv2('D:/growth/digits/data_muzi_all.csv', sep=";")
interp <- growthfd.digits(data=data, colName='X1_2_sur', minCount=9, ncores=4)

resampled <- growthfd.bgs.resample(interp)
smoothed <- growthfd.bgs.smooth(resampled)

age <- seq(10, 18, 0.05)
m<-growthfd.bgs.evalMonotone(smoothed$fd,age,1)
apvs<-growthfd.bgs.apvs(age,m)
ids <- unique(data$id)

age <- seq(min(data$age), max(data$age), 0.05)
values <- growthfd.bgs.evalMonotone(smoothed,age)
vel <- growthfd.bgs.evalMonotone(smoothed,age,1)
acc <- growthfd.bgs.evalMonotone(smoothed,age,2)

growthfd.bgs.plotIndividuals(age, ids, apvs, values, vel, acc, gather, 'plot_X1_2_sur.pdf')
