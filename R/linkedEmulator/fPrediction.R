fPrediction <- function(fModel, fLocsPred, kernel = "mat52") {

  #preprocessing necessary in R so that the matrix concatenates correctly
  if (is.integer(dim(fLocsPred)) == TRUE)
  {
    testTrendMat <- cbind(fLocsPred, 1) #the trend always be this way for f
  } else {
    testTrendMat <- cbind(t(as.matrix(fLocsPred)), 1)
  }
  
  if (kernel == "sq_exp" || kernel == "pow_exp") {
      fPred <- predict.ppgasp(fModel, fLocsPred, testing_trend = testTrendMat, kernel_type = 'pow_exp')
  } else if (kernel == "mat32" || kernel == "matern_3_2") {
        fPred <- predict.ppgasp(fModel, fLocsPred, testing_trend = testTrendMat, kernel_type = 'matern_3_2')
  } else {
        fPred <- predict.ppgasp(fModel, fLocsPred, testing_trend = testTrendMat)
  }
  
  outList <- list("mean" = fPred$mean, "sd" = fPred$sd)
  return(outList) 
}