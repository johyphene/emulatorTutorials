gPrediction <- function(gModel, gLocsPred, trend = "linear", kernel = "mat52") {

	numTest <- nrow(gLocsPred)
	if(trend == 'constant') {
		testTrendMat <- matrix(1,numTest,1)
	} else {
		testTrendMat <- cbind(gLocsPred, matrix(1,numTest,1))
	}
  
	if (kernel == "sq_exp" || kernel == "pow_exp") {
	  gPred <- predict.ppgasp(gModel, gLocsPred, testing_trend = testTrendMat, kernel_type = 'pow_exp')
	} else if (kernel == "mat32" || kernel == "matern_3_2") {
	  gPred <- predict.ppgasp(gModel, gLocsPred, testing_trend = testTrendMat, kernel_type = 'matern_3_2')
	} else {
	  gPred <- predict.ppgasp(gModel, gLocsPred, testing_trend = testTrendMat)
	}

	return(gPred) 
}