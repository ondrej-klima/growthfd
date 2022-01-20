filename <- system.file("extdata", "data.csv", package="growthfd", mustWork=TRUE)
csv <- read.csv(filename)
d <- data.frame('id'=as.factor(csv[,'id']), 'x'=csv[,'age'], 'y'=csv[,'height'])
m <- d$id == 'John'
fit <- growthfd.fit(model.bgs.m, age=d$x[m], height=d$y[m])
data<-growthfd.ApvRegVelocity(model.bgs.m, fit$par)