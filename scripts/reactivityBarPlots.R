library(ggplot2)
library(plyr) 
library(ggthemes)
clusterPlotData = read.csv("clusterPlotData.csv")

predicted = ddply(clusterPlotData, .(variable), function(df) return(df$value > .2))

nCl1 = as.numeric(predicted[1,-1])
nCl2 = as.numeric(predicted[2,-1])
rCl1 = as.numeric(predicted[3,-1])
rCl2 = as.numeric(predicted[4,-1])

nCl1 %*% rCl1 / length(rCl1)
nCl1 %*% rCl2 / length(rCl1)
nCl2 %*% rCl2 / length(rCl1)
nCl2 %*% rCl1 / length(rCl1)

comp1 = subset(clusterPlotData, variable %in% c("rCl1", "nCl1"))
p1 = ggplot(comp1, aes(x=index, y = value, fill = variable)) + geom_bar(stat="identity",position="dodge")
p1 + theme_bw()
ggsave("NestorReactitiyPlot1.eps")

comp2 = subset(clusterPlotData, variable %in% c("rCl2", "nCl2"))
p2 = ggplot(comp2, aes(x=index, y = value, fill = variable)) + geom_bar(stat="identity",position="dodge")
p2 + theme_bw()
ggsave("NestorReactitiyPlot12.eps")

## the scatter plots look really bad...
p3 = qplot(x = subset(comp1, variable == "rCl1")$value, y = subset(comp1, variable == "nCl1")$value)
p3
cor(subset(comp1, variable == "rCl1")$value, subset(comp1, variable == "nCl1")$value)
cor(subset(comp2, variable == "rCl2")$value, subset(comp2, variable == "nCl2")$value)

p4 = qplot(x = subset(comp2, variable == "rCl2")$value, y = subset(comp2, variable == "nCl2")$value)
p4
