cl1 = read.table("cl1.pi")
cl2 = read.table("cl2.pi")
cl3 = read.table("cl3.pi")
cl4 = read.table("cl4.pi")
cl5 = read.table("cl5.pi")
cl6 = read.table("cl6.pi")
cl7 = read.table("cl7.pi")
cl8 = read.table("cl8.pi")
cl9 = read.table("cl9.pi")

react = as.matrix(read.table("bistable.tab")[,-1])
wt = react[1,]
wt = ((1 + tanh((wt - median(wt))^2))/4)
plot_dat = data.frame(index = 1:length(wt), wt = wt)
p = ggplot(data = plot_dat, aes(x = index, y = wt)) + geom_bar(stat = "identity") +
    xlab("Nucleotide Index")+
    ylab("Normalized Reactivity") +
    ggtitle("Experimental Pair Probabilities") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
        )
p

p2_dat = data.frame(index = 1:length(wt), cluster = 1-cl1$V2[1:length(wt)])
p2 = ggplot(data = p2_dat, aes(x = index, y = cluster)) + geom_bar(stat = "identity") +
    xlab("Nucleotide Index")+
    ylab("In-Cluster Probability") +
    ggtitle("Single Cluster Probabilities") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
        )
p2


library(limSolve)
## Use function lsei
## A is a matrix containing clusters
A = cbind(1-cl1$V2, 1-cl2$V2, 1-cl3$V2, 1-cl4$V2, 1-cl5$V2, 1-cl6$V2, 1-cl7$V2, 1-cl8$V2, 1-cl9$V2)[1:length(wt),]
## B is wt (we want to minimize ||Cp - wt||^2)
b = wt
## inequality constraints
## all p's must be positive
g1 = diag(ncol(A))
h1 = rep(0, ncol(A))
# p's must sum to less than 1
g2 = - rep(1, ncol(A))
h2 = -1
G = rbind(g1, g2)
h = c(h1, h2)
sol = lsei(A = A, B = b, G = G, H = h, type = 2)

p3_dat = data.frame(index = 1:length(wt), solution = A%*%sol$X)
p3 = ggplot(data = p3_dat, aes(x = index, y = solution)) + geom_bar(stat = "identity")+
    xlab("Nucleotide Index")+
    ylab("Fit Probability") +
    ggtitle("Best-Fit Probabilities") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
        )
p3

dat_residual = data.frame(index = 1:length(wt), residual = A%*%sol$X - wt)
p4 = ggplot(data = dat_residual, aes(x = index, y = residual)) + geom_bar(stat = "identity") +
    xlab("Nucleotide Index")+
    ylab("Probability Difference") +
    ggtitle("Model Residuals") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
        )
p4


new_react = t(apply(react, 1, function(x) ((1 + tanh((x - median(x))^1))^5/2)))
pca = eigen(cov(new_react))
rotations = pca$vectors[,1:10]
colnames(rotations) = paste("PC", 1:10, sep = "")
## change coordinates to the new ones
projections = as.data.frame(new_react %*% rotations)
clusters = as.data.frame(t(A) %*% rotations)
projections$type = "experiment"
clusters$type = paste("cluster", 1:9, sep = "")

clust_dat = rbind(projections, clusters)
p_clust = ggplot(data=projections, aes( x = PC1, y = PC2)) +
    geom_point() +
    geom_point(data = clusters, color = "red") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
        )
p_clust
