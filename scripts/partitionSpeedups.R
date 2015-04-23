library(ggplot2)
library(ggthemes)

dataNew = read.table("../data/testResultsPFNew.txt", head = TRUE, stringsAsFactors = FALSE)
dataOld = read.table("../data/testResultsPFOld.txt", head = TRUE, stringsAsFactors = FALSE)

dataAll = rbind(dataNew, dataOld)

partitionTimePlot = ggplot(data=dataAll, aes(x=n, y=time, color=type)) + geom_point()