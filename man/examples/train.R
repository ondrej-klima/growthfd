data <- read.table("D:/Growth/playground/data/bgs_07.txt", header=TRUE, sep="\t", na.strings="NA", dec=".")
data <- data[!(data$id %in% growthfd.bgs.dropoutsIds.Height()) & data$sex == "1",]
gather <- growthfd.bgs.gather(data)
interp <- growthfd.bgs.interpolateNAs(gather)
resampled <- growthfd.bgs.resample(interp)
smoothed <- growthfd.bgs.smooth(resampled)

age <- seq(10, 18, 0.05)
m<-growthfd.bgs.evalMonotone(smoothed,age)
apvs<-growthfd.bgs.apvs(age,m)

age <- seq(0, 18, 0.05)
ids <- unique(data$id)
values <- growthfd.bgs.evalMonotone(smoothed,age)
vel <- growthfd.bgs.evalMonotone(smoothed,age,1)
acc <- growthfd.bgs.evalMonotone(smoothed,age,2)
growthfd.bgs.plotIndividuals(age, ids, apvs, values, vel, acc, gather)
growthfd.bgs.plotAll(age,acc, ylimit = c(-25, 25))