filename <- system.file("extdata", "data.csv", package="growthfd", mustWork=TRUE)
csv <- read.csv(filename)
d <- data.frame('id'=as.factor(csv[,'id']), 'x'=csv[,'age'], 'y'=csv[,'height'])
growthfd(data=d, x=x, y=y, id=id, model=model.bgs.m)