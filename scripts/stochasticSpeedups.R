library(ggplot2)

data1 = read.table("../data/test_results1.txt", head = TRUE, stringsAsFactors = FALSE)
data2 = read.table("../data/test_results2.txt", head = TRUE, stringsAsFactors = FALSE)

data1$type[which(data1$type == "stochastic")] = "oldstochastic"

timingData = rbind(data1, data2)
timingData = subset(timingData, type != "partition")

p = ggplot(data = timingData, aes(x = n, y = time, color = type)) + geom_point()
