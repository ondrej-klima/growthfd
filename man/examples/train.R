data <- read.table("D:/Growth/playground/data/bgs_07.txt", header=TRUE, sep="\t", na.strings="NA", dec=".")
data <- data[!(data$id %in% growthfd.bgs.dropoutsIds.Height()) & data$sex == "1",]
gather <- growthfd.bgs.gather(data)
interp <- growthfd.bgs.interpolateNAs(gather)
resampled <- growthfd.bgs.resample(interp)
smoothed <- growthfd.bgs.smooth(resampled)

age <- seq(10, 18, 0.05)
m<-growthfd.bgs.evalMonotone(smoothed,age,1)
apvs<-growthfd.bgs.apvs(age,m)

age <- seq(0, 18, 0.05)
ids <- unique(data$id)
values <- growthfd.bgs.evalMonotone(smoothed,age)
vel <- growthfd.bgs.evalMonotone(smoothed,age,1)
acc <- growthfd.bgs.evalMonotone(smoothed,age,2)
growthfd.bgs.plotIndividuals(age, ids, apvs, values, vel, acc, gather)
growthfd.bgs.plotAll(age,acc, ylimit = c(-25, 25))

smoothed2 <- growthfd.bgs.smooth(resampled, F)
m2<-growthfd.bgs.eval(smoothed2$fd,age)
values2 <- growthfd.bgs.eval(smoothed2$fd,age)
vel2 <- growthfd.bgs.eval(smoothed2$fd,age,1)
acc2 <- growthfd.bgs.eval(smoothed2$fd,age,2)
growthfd.bgs.plotIndividuals(age, ids, apvs, values2, vel2, acc2, gather, filename = "plots2.pdf")

resampled3 <- values
dim(resampled3) <- c(prod(dim(values)), 1)
resampled3 <- cbind(rep(ids, each=dim(values)[1]), rep(age, dim(values)[2]), resampled3[,1])
colnames(resampled3) <- colnames(resampled)
smoothed3 <- growthfd.bgs.smooth(resampled3, F)
m3<-growthfd.bgs.eval(smoothed3$fd,age)
values3 <- growthfd.bgs.eval(smoothed3$fd,age)
vel3 <- growthfd.bgs.eval(smoothed3$fd,age,1)
acc3 <- growthfd.bgs.eval(smoothed3$fd,age,2)
growthfd.bgs.plotIndividuals(age, ids, apvs, values3, vel3, acc3, gather, filename = "plots3.pdf")

tw <- growthfd.bgs.registerCurvesToApvs(smoothed, apvs)
reged<-fda::register.newfd(smoothed$yhatfd, tw)
fda::register.newfd(accelfdUN, regListLM$warpfd)

reged.bak<-fda::register.newfd(smoothed$yhatfd, tw)
tw.bak <- tw

amp<-fda::register.newfd(smoothed3$fd,tw)
itw <- growthfd.bgs.invertTw(age, tw)

model <- growthfd.bgs.model(amp, itw$fd)