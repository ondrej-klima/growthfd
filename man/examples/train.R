data <- read.table("D:/Growth/playground/data/bgs_07.txt", header=TRUE, sep="\t", na.strings="NA", dec=".")
data <- data[!(data$id %in% growthfd.bgs.dropoutsIds.Height()) & data$sex == "1",]
gather <- growthfd.bgs.gather(data)
interp <- growthfd.bgs.interpolateNAs(gather)
resampled <- growthfd.bgs.resample(interp)
smoothed <- growthfd.bgs.smooth(resampled)
