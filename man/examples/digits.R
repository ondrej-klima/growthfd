data <- read.csv2('D:/growth/digits/data_muzi_upravene03.csv', sep=";")
#data <- read.csv2('D:/growth/digits/data_zeny_all.csv', sep=";")
colNames <- c('mt1','mt2','mt3','mt4','mt5','pp1','pp2','pp3','pp4','pp5','mp2',
              'mp3','mp4','mp5','dp1','dp2','dp3','dp4','dp5','f1','f2','f3',
              'f4','f5','pw15','dw15','pw25','dw25','lh', 'height')
colNames <- c('height')

#data <- read.csv2('D:/growth/digits/data_muzi_all.csv', sep=";")
colNames<-c('mt2_A', 'mt2_B')
colNames<-c('mt2_B')

dropIds1 <- c(1002, 1007, 1023, 1024, 1027, 1038, 1044, 1047, 1053, 1062, 1076,
              1100, 1111, 1130, 1133, 1137, 1138, 1139, 1143, 1160, 1175, 1178, 
              1190, 1198, 1205, 1225, 1230, 1232, 1242, 1250, 1252, 1277, 1290,
              1299, 1327, 1329, 1330, 1344, 1357, 1360, 1425)

dropIds2 <- c(1009, 1054, 1055, 1069, 1071, 1161, 1169, 1221, 1343, 1346, 1386,
              1398)

for(colName in colNames) {
  #res <- growthfd.digits(data=data, colName=colName, minCount=9, ncores=4)
  res <- growthfd.digits.spline(data=data, colName=colName, minCount=9)
  interp <- res[[1]]
  interp$value <- interp$valuei
  
  gather <- as.data.frame(data[,c('ind', 'age', colName)])
  colnames(gather) <- c('id', 'age', 'value')
  
  resampled <- growthfd.bgs.resample(interp)
  #smoothed <- growthfd.bgs.smooth(resampled, F)
  smoothed <- growthfd.bgs.smooth(resampled, norder=8)
  
  age <- seq(10, 18, 0.05)
  #m<-growthfd.bgs.eval(smoothed$fd,age,1)
  m<-growthfd.bgs.evalMonotone(smoothed,age,1)
  apvs<-growthfd.bgs.apvs(age,m)
  ids <- unique(interp$id)
  
  age<-seq(min(interp$age), max(interp$age), 0.05)
  
  values <- growthfd.bgs.evalMonotone(smoothed,age)
  vel <- growthfd.bgs.evalMonotone(smoothed,age,1)
  acc <- growthfd.bgs.evalMonotone(smoothed,age,2)
  
  #values <- growthfd.bgs.eval(smoothed$fd, age)
  #vel <- growthfd.bgs.eval(smoothed$fd,age,1)
  #acc <- growthfd.bgs.eval(smoothed$fd,age,2)
  
  growthfd.bgs.plotIndividuals(age, ids, apvs, values, vel, acc, gather, sprintf('plots/mt2_b/%s_9_2_n_mon.pdf', colName))
  growthfd.bgs.plotIndividualsPoints(unique(res[[2]]), gather,  sprintf('plots/mt2_b/%s_pts_9_2_n_mon.pdf', colName))
}

# Fit ordinary splines to monotone splines
resampled.new <- values
dim(resampled.new) <- c(prod(dim(values)), 1)
resampled.new <- cbind(rep(ids, each=dim(values)[1]), rep(age, dim(values)[2]), resampled.new[,1])
colnames(resampled.new) <- colnames(resampled)
smoothed.new <- growthfd.bgs.smooth(resampled.new, F)

tw <- growthfd.bgs.registerCurvesToApvs(smoothed, apvs)
amp<-fda::register.newfd(smoothed.new$fd,tw)
itw <- growthfd.bgs.invertTw(age, tw)

model <- growthfd.bgs.model(amp, itw$fd)
growthfd.plotwarps(model, age=seq(7.5, 20, 0.05))

data.df <- as.data.frame(data)
data.df$id <- as.factor(data.df$ind)
#data.df$value <- data.df$mt2_A
data.df$value <- data.df$mt2_B

#r2 <- growthfd(data=data.df[data.df$id %in% c(1001,1002,1003,1004),], x=age, y=mt2_A, model = model, id=id, parallel = T)
r <- growthfd(data=data.df, x=age, y=mt2_B, model = model, id=id, parallel = T)
growthfd.bgs.plotIndividuals(age = r$sampling, ids = r$ids, apvs = r$milestones[,1], values=t(r$stature), vel = t(r$velocity), acc=t(r$acceleration), data = data.df, filename='plot_mt2_B_fit.pdf')

cor <- cbind(as.numeric(as.character(r$fitted[,'id'])), r$fitted[,'age']-r$fitted[,'apv'], r$fitted[,'residuum'])
colnames(cor) <- c('id', 'age', 'res')
write.csv2(cor, file='residuals.on.apv.mt2B.csv')

x<-data$age[data$ind==1017]
y<-data$mt2_B[data$ind==1017]
sp<-smooth.spline(x,y,df=length(x)-1)
sq<-seq(10,20,0.05)
plot(sq, predict(sp,sq)$y)
points(x,y,pch=4)


