data <- read.csv2('D:/growth/digits/data_muzi_all.csv', sep=";")
res <- growthfd.digits(data=data, colName='X1_2_sur', minCount=1, ncores=4)
interp <- res[[1]]
interp$value <- interp$valuei

gather <- as.data.frame(data[,c('ind', 'age', 'X1_2_sur')])
colnames(gather) <- c('id', 'age', 'value')

resampled <- growthfd.bgs.resample(interp)
smoothed <- growthfd.bgs.smooth(resampled, F)

age <- seq(10, 18, 0.05)
m<-growthfd.bgs.eval(smoothed$fd,age,1)
apvs<-growthfd.bgs.apvs(age,m)
ids <- unique(interp$id)

age<-seq(min(interp$age), max(interp$age), 0.05)
values <- growthfd.bgs.eval(smoothed$fd,age)
vel <- growthfd.bgs.eval(smoothed$fd,age,1)
acc <- growthfd.bgs.eval(smoothed$fd,age,2)
growthfd.bgs.plotIndividuals(age, ids, apvs, values, vel, acc, gather, 'plot_X1_2_sur.pdf')
growthfd.bgs.plotIndividualsPoints(unique(res[[2]]), gather, 'plot_X1_2_sur_pts.pdf')