linkedEmulator <- function(fModels, gModel, fLocsPred, zLocsPred, trend = "linear", kernel = "mat52", eta = 0.000000000001) {
  #load in correct .cpp file for kernel specified
  library(Rcpp)
  setwd("C:/Users/jl541227/Documents/ResearchPapers/Richards_Joey2/Richards_Joey/")
  
  if (kernel == "exp") {
    sourceCpp('linkedEmulatorExponential.cpp') # file with C++ versions of the linked emulator functions for exponential
  } else if (kernel == "sq_exp" || kernel == "pow_exp") {
    sourceCpp('linkedEmulatorSquaredExponential.cpp') # file with C++ versions of the linked emulator functions for squared exponential
  } else if (kernel == "mat32") {
    sourceCpp('linkedEmulatorMatern32.cpp') # file with C++ versions of the linked emulator functions for Matern 3/2
  } else {
    sourceCpp('linkedEmulator.cpp') # file with C++ versions of the linked emulator functions for Matern 5/2
  }
  
  set.seed(9)
  
  N <- nrow(as.matrix(fLocsPred))
  numFInputs <- dim(fModels[[1]]@input)[2]
  numFOutputs <- length(fModels)
  ygpVals <- matrix(0, nrow = N, ncol = numFOutputs)
  sdVals <- matrix(0, nrow = N, ncol = numFOutputs)
  gammas <- matrix(0, nrow = numFOutputs, ncol = numFInputs)
  fGamma <- rep(1, times=numFInputs)
  
  #first layer PPE - take in f^D inputs/outputs, take in f^* inputs
  #a single emulator is constructed for each output variable if we have more than
  #one. Predictions are more accurate this way.
  for (jj in 1:numFOutputs) {
    fVals = fPrediction(fModels[[jj]], fLocsPred, kernel)
    
    ygp = fVals$mean
    v = fVals$sd
    ygpVals[,jj] <- ygp
    sdVals[,jj] <- v
    gammas[jj,] <- 1/fModels[[jj]]@beta_hat
    for (jjj in 1:numFInputs) {
      fGamma[jjj] <- fGamma[jjj]*gammas[jj,jjj]
    }
  }

  fPred <- list(mean = ygpVals, sd = sdVals)
  
  gLocs <- gModel@input
  
  gLocsPred <- cbind(fPred$mean, zLocsPred)

  
  #train inputs are called gLocs, train outputs are called gOutput, test inputs are called fLocsPred
  gPred = gPrediction(gModel, gLocsPred, trend, kernel)
  
  muVec = as.matrix(fPred$mean)
  sigVec = as.matrix(fPred$sd)
  
  #number of inputs from f is determined by number of columns in muVec - PPE mean prediction
  #number of external inputs for g determined by subtracting the f input number from the total in g
  d <- ncol(muVec)
  p <- ncol(gLocs) - d
  wT <- gLocs[,1:d]
  zT <- gLocs[,(d+1):(d+p)]
  m <- dim(gLocs)[1]
  gSigVec <- gModel@sigma2_hat
  gGamma <- 1/gModel@beta_hat
  gOutput <- gModel@output

  #calculation correlations between training data points
  c_k <- array(0, dim = c(m, m, (d+p)))
  for (i in 1:m) {
    for (j in 1:m) {
      for (k in 1:(d+p)) {
        c_k[i,j,k] <- cKSingle(gLocs[i,k], gLocs[j,k], gGamma[k])
      }
    }
  }
  for (k in 1:(d+p)) {
    if (k == 1)
      R <- c_k[,,k]    
    else 
      R <- R * c_k[,,k]
  }
  
  #adding term here to help with regularizing before inverting
  invR <- solve(R + diag(eta, nrow(R), ncol(R)))

  g = ncol(as.matrix(gOutput)) #number of outputs we are interested in

  numTest <- nrow(gLocsPred)
  
  #function gets split into 4 options from here:
  #one output with either a linear or constant trend 
  #OR
  #multiple outputs with either a linear or constant trend
  #each of these cases has its own function call which contains the linkedSampling function
  if (g == 1) { #if we only have one output of interest
    if (trend == "constant") { #call the constant trend, single output case
      samps <- linkedConstantSingle(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta)
    } else { #all other entries will default to the linear case - so call the linear trend, single output case
      samps <- linkedLinearSingle(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta)
    }
  } else { #g > 1
      if (trend == "constant") { #call the constant trend, multiple output case
      samps <- linkedConstantMultiple(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta)
    } else { #all other entries will default to the linear case - so call the linear trend, multiple output case
      samps <- linkedLinearMultiple(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta)
    }
  }

  outList <- list("means" = gPred$mean, "vars" = samps$A, "samples" = samps$samples)
  return(outList)
}