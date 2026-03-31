uniqueRows = function(m){
  # returns row numbers of a set of unique rows (not a matrix)
  m = as.matrix(m)
  nr = dim(m)[1]
  nc = dim(m)[2]
  m2 = data.matrix(m)
  unifRandomVec = runif(nc)
  v = matrix(unifRandomVec,nc,1)
  r = m %*% v
  idx = sort(r,index.return = TRUE)
  idx = idx[[2]]
  rowidx = idx[1]
  for (i in 2:nr) {
    if (abs(r[idx[i]] - r[idx[i - 1]]) > .1e-6) {
      rowidx = rbind(rowidx,idx[i])
    }
  }
  rowidx = sort( rowidx )
  return(rowidx)
}


calcPowersForRepPatterns = function(U,n_T,mu) {
  t = nrow(U)
  nReps = n_T - t
  QU = U %*% solve(t(U) %*% U) %*% t(U)
  eig = eigen(QU)
  cols = which(abs(eig$values) < .1e-10)
  U0 = eig$vectors[ ,cols]
  gridList = as.list(1:nReps)
  for (i in 1:nReps) {
    gridList[[i]] = 1:t
  }
  grid = expand.grid(gridList)
  cases = matrix(0,nrow(grid),t)
  for (i in 1:nrow(grid)) {
    for (j in 1:t) {
      cases[i,j] = sum(grid[i,] == j)
    }
  }
  allCases = cases
  keepers = uniqueRows(cases)
  cases = cases[keepers, ]
  nCases = nrow(cases)
  traceMat = matrix(0,nCases,1)
  for (i in 1:nCases) {
    W = diag(1,t,t) + diag(cases[i, ],t,t)
    Winv = solve(W)
    traceMat[i] = sum(diag(solve(t(U0) %*% Winv %*% U0)))
  }
  cases = cbind(cases,traceMat)
  casesPrior = matrix(0,2,ncol(cases))
  Uproj = U %*% solve(t(U) %*% U) %*% t(U)
  tau = mu - Uproj %*% mu
  casesPrior[1,1:t] = tau
  casesPrior[2,1:t] = diag(Uproj)
  return(list(cases,casesPrior))
}


U = as.matrix(expand.grid(x = c(-1,0,1),y = c(-1,0,1)))
U = cbind(matrix(1,t,1),U[,1],U[,2])
for (analNum in 1:5) {
  t = nrow(U)
  #
  # Case 1-5
  #
  if (is.element(analNum,c(1,4,5))) {
    mu = 2 + U[,2] + U[,2]^2 + U[,3] + U[,3]^2 + rnorm(nrow(U),0,.2)
    # Outlier for analysis 5
    if (analNum == 5) {
      for (iTrt in 1:t) {
        if ((U[iTrt,2] + U[iTrt,3]) == 2) {
          mu[iTrt] = mu[iTrt] + 2
        }
      }
    }
  }
  if (analNum == 2) {
    mu = 2 + U[,2] + U[,2]^2 + U[,3] + U[,3]^2 + rnorm(nrow(U),0,.6)
  }
  #
  # Case 3
  #
  if (analNum == 3) {
    mu = 2 + U[,2] + U[,2]^2 + U[,3] + U[,3]^2 + rnorm(nrow(U),0, .4*(U[,2]^2 + U[,3]^2))
  }
  nCols = dim(U)[2]
  if (analNum < 4) {
    n_T = t + 3
  } else {
    n_T = t + 5
  }
  casesOut = calcPowersForRepPatterns(U,n_T,mu)
  results = casesOut[[1]]
  results = data.frame(results)
  results <- results[order(results[, ncol(results)],decreasing = TRUE), ]
  colnames(results) = c(paste0("X",i = 1:t),"lamda")
  tauAndH = data.frame(casesOut[[2]])
  colnames(tauAndH) = c(paste0("X",i = 1:t),"lamda")
  results = rbind(tauAndH,results)
  nRows = nrow(results)
  results = rbind(results[1:27, ],results[(nRows - 24):nRows, ])
  if (analNum == 1) {
    fName = paste0("Delta1.xlsx")
  }
  if (analNum == 2) {
    fName = paste0("Delta2.xlsx")
  }
  if (analNum == 3) {
    fName = paste0("Delta3.xlsx")
  }
  if (analNum == 4) {
    fName = paste0("Delta4.xlsx")
  }
  if (analNum == 5) {
    fName = paste0("Delta5.xlsx")
  }
  write.xlsx(results,fName)
}



