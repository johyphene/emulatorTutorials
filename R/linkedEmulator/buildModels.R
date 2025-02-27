buildModels <- function(fLocs, fOutput, gLocs, gOutput, trend, kernel) {

	#we build the f model first
	m <- nrow(fLocs)
	fTrendMat <- cbind(fLocs, matrix(1,m,1)) #the trend always be this way for f
	
	M <- nrow(gLocs)
	if(trend == 'constant') {
		gTrendMat <- matrix(1,M,1)
	} else {
    gTrendMat <- cbind(gLocs, matrix(1,M,1))
	}
	
	fModels <- list()
	
	#an individual f emulator is created for each output component to improve prediction accuracy
	
	#determining size of output and then using the correct optimization based on size
	fSwitchVal = dim(as.matrix(fOutput))[1] * dim(as.matrix(fOutput))[2]
	for	(i in 1:dim(as.matrix(fOutput))[2]) {
		if (fSwitchVal < 10000) {
			if (kernel == "sq_exp" || kernel == "pow_exp") {
				fModel <- ppgasp(design=fLocs, response=fOutput[,i], trend=fTrendMat, kernel_type = 'pow_exp') #ppgasp for f layer
			} else if (kernel == "mat32") {
				fModel <- ppgasp(design=fLocs, response=fOutput[,i], trend=fTrendMat, kernel_type = 'matern_3_2') #ppgasp for f layer
			} else {
				fModel <- ppgasp(design=fLocs, response=fOutput[,i], trend=fTrendMat) #ppgasp for f layer 
			}
		} else {
			if (kernel == "sq_exp" || kernel == "pow_exp") {
				fModel <- ppgasp(design=fLocs, response=fOutput[,i], trend=fTrendMat, kernel_type = 'pow_exp', optimization = "nelder-mead") #ppgasp for f layer
			} else if (kernel == "mat32") {
				fModel <- ppgasp(design=fLocs, response=fOutput[,i], trend=fTrendMat, kernel_type = 'matern_3_2', optimization = "nelder-mead") #ppgasp for f layer
			} else {
				fModel <- ppgasp(design=fLocs, response=fOutput[,i], trend=fTrendMat, optimization = "nelder-mead") #ppgasp for f layer
			}
		}
		fModels[[i]] <- fModel
	}
	
	#determining size of output and then using the correct optimization based on size
	gSwitchVal = dim(as.matrix(gOutput))[1] * dim(as.matrix(gOutput))[2]
	if (gSwitchVal < 10000) {
		if (kernel == "sq_exp" || kernel == "pow_exp") {
			gModel <- ppgasp(design=gLocs, response=gOutput, trend=gTrendMat, kernel_type = 'pow_exp') #ppgasp for g layer
		} else if (kernel == "mat32") {
			gModel <- ppgasp(design=gLocs, response=gOutput, trend=gTrendMat, kernel_type = 'matern_3_2') #ppgasp for g layer
		} else {
			gModel <- ppgasp(design=gLocs, response=gOutput, trend=gTrendMat) #ppgasp for g layer
		}
	} else {
		if (kernel == "sq_exp" || kernel == "pow_exp") {
			gModel <- ppgasp(design=gLocs, response=gOutput, trend=gTrendMat, kernel_type = 'pow_exp', optimization = "nelder-mead") #ppgasp for g layer
		} else if (kernel == "mat32") {
			gModel <- ppgasp(design=gLocs, response=gOutput, trend=gTrendMat, kernel_type = 'matern_3_2', optimization = "nelder-mead") #ppgasp for g layer
		} else {
			gModel <- ppgasp(design=gLocs, response=gOutput, trend=gTrendMat, optimization = "nelder-mead") #ppgasp for g layer
		}
	}
	outlist <- list("fModels" = fModels, "gModel" = gModel)
	return(outlist)
}