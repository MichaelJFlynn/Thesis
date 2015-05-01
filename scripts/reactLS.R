

# probably should use the optim function

## Computes the cost function for the parameters
## pmus should be matrix of nest probabilities
## pnus should be the in-nest base probabilities
## for k clusters, n bases, and m experiments:
## pmus is m x k
## pis is k x n
## reacts is m x n
## In par, the first m*k parameters go in pmus
## the next n*k go to pis
cost_function <- function(par, reacts, k) {
  defect = 0;
  ## constraints: all rows of pmus must sum to at most 1
  m = nrow(reacts)
  n = ncol(reacts)
  rowsums = rowSums(matrix(par[1:(m*k)], ncol = k))
  defect = defect + sum(rowsums[which(rowsums > 1)]-1) 
  stopifnot(length(par) == k*(n+m))
  pmus = matrix(par[1:(k*nrow(reacts))],
    nrow = nrow(reacts),
    ncol = k)
  pnuis = matrix(par[(k*nrow(reacts)+1):length(par)],
    nrow = k,
    ncol = ncol(reacts))
  pis = pmus %*% pnuis;
  costs = pis - reacts;
  return(sum(apply(costs, 1, function(x) x %*% x))*exp(defect));
}

## since matrix is filled up column by column, we want to take it
## apart similarly
mat_to_par = function(pmus, pnuis) {
  return(c(pmus, pnuis))
}

par_to_mat = function(par, k, n, m) {
  stopifnot(length(par) == k*(n+m))
  pmus = matrix(par[1:(k*m)],
    nrow = m,
    ncol = k)
  pnuis = matrix(par[(k*m+1):(k*(n+m))],
    nrow = k,
    ncol = n)
  return(list(pmus = pmus, pnuis = pnuis))
}

k = 5
m = 20
n = 60
pmus = matrix(rnorm(k*m), nrow = m, ncol = k)
pnuis = matrix(rnorm(k*n), nrow = k, ncol = n)
reacts = matrix(rnorm(m*n), nrow = m, ncol = n)
cost_function(mat_to_par(pmus, pnuis), reacts, k)
tests = par_to_mat(mat_to_par(pmus, pnuis), k=k, n=n, m=m)
## assert test cases:
all(tests$pmus == pmus)
all(tests$pnuis == pnuis)

## now get reactivity data
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
n01pairs = read.table("../n01.pairs", skip = 1)
n02pairs = read.table("../n02.pairs", skip = 1)

## transpose because apply default returns by column
normalized = t(apply(react, 1, function(x) 1 - (1+tanh(x - median(x)))^2/4))

## now fit a model using the clusters as starting point
## also cut off some data because we're missing the last 20 BPs
pnuis = rbind(cl1$V2, cl2$V2, cl3$V2, cl4$V2, cl5$V2, cl6$V2, cl7$V2, cl8$V2, cl9$V2)[,1:ncol(react)]
pnuis = pnuis[1:2,]
k = nrow(pnuis)
n = ncol(normalized)
m = nrow(normalized)
pmus = matrix(1/k, ncol = k, nrow = m)

## optimize!
sol = optim(par = mat_to_par(pmus, pnuis),
  fn = cost_function,
  method = "L-BFGS-B",
  lower = 0,
  upper = 1,
  k = k,
  react = normalized,
  control = list(maxit = 10000))
sol_mats = par_to_mat(sol$par, k=k, n=n, m=m)

## Plots to analyze, looks close
library(ggplot2)
qplot(x = 1:n, y = c(sol_mats$pmus[1,] %*% sol_mats$pnuis), geom="bar", stat = "identity")
qplot(x = 1:n, y = normalized[1,], geom = "bar", stat = "identity")
qplot(x = 1:n, y = c(sol_mats$pmus[1,] %*% sol_mats$pnuis) - normalized[1,], geom="bar", stat = "identity")

## now, compare to nestor
apply(sol_mats$pnuis, 1, function(x) apply(pnuis, 1, function(y) sum((x-y)^2))) -> diffs
apply(sol_mats$pnuis, 1, function(x) order(apply(pnuis, 1, function(y) (x-y)%*%(x-y))))
result_clusts = sol_mats$pmus %*% sol_mats$pnuis
apply(normalized-result_clusts, 1, function(x) sum(x^2))

library(limSolve)
A = t(pnuis)
G = rbind(diag(ncol(A)), rep(-1, ncol(A)))
H = c(rep(0,ncol(A)), -1)
## check Nestor cluster solutn
apply(normalized, 1, function(wt) {
  b = wt
  sol = lsei(A = A, B = b, G = G, H = H, type = 2)
  return(sum((A%*%sol$X - b)^2))
})

wt = normalized[1,]
first = result_clusts[1,]
## sum of squares pretty close
(wt-first)%*%(wt-first)

## plot all the k clusters
fit_clusters = as.data.frame(t(sol_mats$pnuis))
colnames(fit_clusters) = paste("rCl", 1:nrow(sol_mats$pnuis), sep = "")
fit_clusters$index = 1:nrow(fit_clusters)
library(reshape)
melted = melt(fit_clusters, id.vars=  "index")
p = ggplot(melted, aes(x = index, y = value, color = variable)) + geom_point() + geom_line()
p

p2 = ggplot(melted, aes(x = index, y= value, fill = variable)) + geom_bar(stat = "identity", position = "dodge")
p2

## now do the same with first 2 nestor clusters?
fit_clusters = as.data.frame(t(pnuis[c(1,2),]))
colnames(fit_clusters) = paste("nCl", 1:nrow(sol_mats$pnuis), sep = "")
fit_clusters$index = 1:nrow(fit_clusters)
library(reshape)
melted = rbind(melted, melt(fit_clusters, id.vars=  "index"))
p = ggplot(melted, aes(x = index, y = value, color = variable)) + geom_point() + geom_line()
p

apply(normalized, 2, mean)

each_variable = cast(melted)

ggplot(subset(melted, variable %in% c("nCl1", "rCl1")), aes(x = index, y= value, fill = variable)) + geom_bar(stat = "identity", position = "dodge")
ggsave("NvsRCluster1.png")
ggplot(subset(melted, variable %in% c("nCl2", "rCl2")), aes(x = index, y= value, fill = variable)) + geom_bar(stat = "identity", position = "dodge")
ggsave("NvsRCluster2.png")

mutateChange = data.frame(sol_mats$pmus)
colnames(mutateChange) = c("CL1", "CL2")
mutateChange$index = 1:nrow(mutateChange)
meltMC = melt(mutateChange, id.vars = "index")
pmC = ggplot(meltMC, aes(x = index, y = value, color = variable)) + geom_line()

## analysis on the two cluster system
## with n01pairs and n02pairs

mnmData = read.table("test_results.txt", stringsAsFactors = F)

library(stringr)
mnmData[,1] = as.numeric(str_extract(mnmData[,1], "[0-9]+"))
mnmData[(is.na(mnmData[,1])), 1] = 0
mnmData[,1] = mnmData[,1] + 1
mnmData[,2] = mnmData[,2]/1000
mnmData[,3] = mnmData[,3]/1000
colnames(mnmData) = c("index", "CL1", "CL2")
meltMC2 = melt(mnmData, id.vars= "index")

meltMC$type = "reactivity fit"
meltMC2$type = "Nestor"

pmC2 = ggplot(rbind(meltMC, meltMC2), aes(x = index, y = value, linetype = type, color = variable)) + geom_line()
