linkedSampling <- function(fLocsPred, gLocsPred, muLs, sigma2Ls, fGamma, gGamma, N, gPred){
  
  A=matrix(1, nrow = N, ncol = N)
  numFCols = ncol(fLocsPred)
  if (is.null(numFCols) == TRUE) {
    numFCols = length(fLocsPred)
    fLocsPred = matrix(fLocsPred, nrow = 1)
  }
  numGCols = ncol(gLocsPred)
  
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) {
        A[i,j] = sigma2Ls[i]
      }
      else if (i != j) {
        rho_f = 1
        for (k in 1:numFCols)  {
          if (nrow(fLocsPred) > 1) {
            rho_f = rho_f * cKSingle(fLocsPred[i,k],fLocsPred[j,k],fGamma[k]) 
            rho_f = rho_f * cKSingle(fLocsPred[k],fLocsPred[k],fGamma[k]) 
          }
        }
        tempSig = sqrt(sigma2Ls[i] - 2*rho_f*sqrt(sigma2Ls[i])*sqrt(sigma2Ls[j]) + sigma2Ls[j])
        rho_eta = 1
        for (k in 1:numGCols)  {
          rho_eta = rho_eta * xiSingle(muLs[i], tempSig, gGamma[k], muLs[j])
        }
        A[i,j] = sqrt(sigma2Ls[i]) * sqrt(sigma2Ls[j]) * rho_eta * rho_f
      }
    }
  }
  L = chol(A + eye(N)*0.0000000001) #regularizing the matrix so it is positive definite
  L = t(L)
  samples <- matrix(0, nrow = N, ncol = 100)
  for (i in 1:100) {
    u = rnorm(N)
    ysamp = gPred + L%*%as.matrix(u)
    samples[,i] <- ysamp
  }
  outList <- list("A" = A, "samples" = samples) 
  return(outList)
  
}