reacts = read.table("bistable.tab")[,-1]
cl1 = read.table("compat.pi")
cl2 = read.table("noncompat.pi")


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
fullPartitionProbabilities = .00876194*cl1$V2 + .991242*cl2$V2
partitionProbsRankSorted = data.frame(rank = 1:nrow(cl1), prob = fullPartitionProbabilities[order(fullPartitionProbabilities, decreasing = T)])
qplot(data = partitionProbsRankSorted, x = rank, y = prob)

map_reacts = function(reacts) {
  
  
}


model = nls(dist~ 1/exp((index- A)/ B))


## now fit a model using the clusters as starting point
## also cut off some data because we're missing the last 20 BPs
pnuis = rbind(cl1$V2, cl2$V2)[,1:ncol(react)]
k = nrow(pnuis)
n = ncol(normalized)
m = nrow(normalized)
pmus = matrix(1/k, ncol = k, nrow = m)

