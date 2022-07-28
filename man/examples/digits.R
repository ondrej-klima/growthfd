data <- read.csv2('D:/growth/digits/data_muzi_all.csv', sep=";")
result <- growthfd.digits(data=data, colName='X1_2_sur', minCount=9, ncores=4)