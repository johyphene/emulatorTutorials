linkedConstantSingle <- function(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta) {
	yT <- gOutput
	HzT <- matrix(1,m,1)
	hTilde <- matrix(1,m,1)
	hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
	thetaHat <- rep(0, times = ncol(as.matrix(wT)))
	thetaHat <- as.matrix(thetaHat)
	betaHat <- as.matrix(hats)
	sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde)
	B <- matrix(0,nrow = d, ncol = m)
	Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
	C <- solve(t(hTilde) %*% invR %*% hTilde)
	muLs <- matrix(c(0), nrow = numTest, ncol = g)
	sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
	
	#calculate all progress locations
	progLocs <- ceil((1:10*(0.1))*(g*numTest))
	
	for (i in 1:g) {
		A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)
		for(k in 1:numTest){
			currLoc <- ((i-1)*numTest) + k
			if(sum(currLoc == progLocs) == 1) { #print statement if 10/20/.../100% through calculations
				print(paste0(floor((currLoc/(g*numTest))*100), "% done"))
			}
			Sys.sleep(0.001) #can remove this when package is ready to go live probably
			flush.console() #can remove this when package is ready to go live probably
			hZ <- 1

			iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			K <- cbind(t(B), iVec%*%t(hZ))
			J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			tempOmega <- sigVec[k,]^2
			Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))
			P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
			G <- cbind(t(muVec[k,]), t(hZ))
			G <- t(G)
			trQJ = 0
			for (bb in 1:ncol(Q)) {
				trQJ = trQJ + (Q[bb,] %*% J[,bb])
			}
			mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
			sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + gSigVec[i]*(1 + eta + sigInvMat + trQJ - 2*tr(C * t(hTilde) %*% invR %*% iVec))
			muLs[k, i] <- mu_L
			sigma2Ls[k, i] <- sigma2_L
			if(sigma2_L < 0) { 
				if(sigma2_L < -(10^-7)) {
					if(sigma2_L < -(10^-7)) {
						print(paste0("Negative variance at index ", k))
						print(paste0("Value of ", sigma2_L))
					}
				}
				sigma2Ls[k,i] = abs(sigma2_L)
			}
		}
	}
	samps = linkedSampling(fLocsPred, gLocsPred, muLs, sigma2Ls, fGamma, gGamma, numTest, gPred$mean)
	return(samps)
}
###
linkedLinearSingle <- function(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta) {
	yT <- gOutput
	HzT <- cbind(zT,matrix(1,m,1))
	hTilde <- cbind(wT, HzT)
	hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% yT
	thetaHat <- as.matrix(hats[1:d,])
	betaHat <- as.matrix(hats[(d+1):(d+p+1),]) 
	sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde)
	B <- matrix(0,nrow = d, ncol = m)
	Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
	C <- solve(t(hTilde) %*% invR %*% hTilde)	
	muLs <- matrix(c(0), nrow = numTest, ncol = g)
	sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
	
	#calculate all progress locations
	progLocs <- ceil((1:10*(0.1))*(g*numTest))
	
	for (i in 1:g) { #should be able to cut this given that g == 1 for this case
		A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)
		for(k in 1:numTest){
			currLoc <- ((i-1)*numTest) + k
			if(sum(currLoc == progLocs) == 1) { #print statement if 10/20/.../100% through calculations
				print(paste0(floor((currLoc/(g*numTest))*100), "% done"))
			}
			Sys.sleep(0.001) #can remove this when package is ready to go live probably
			flush.console() #can remove this when package is ready to go live probably
			hZ <- c(gLocsPred[k,(d+1):(d+p)],1)

			iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			K <- cbind(t(B), iVec%*%t(hZ))
			J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			tempOmega <- sigVec[k,]^2
			if (length(tempOmega) == 1) {
				Omega <- as.matrix(tempOmega)
			} else { 
				Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))
			}
			P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
			G <- cbind(t(muVec[k,]), t(hZ))
			G <- t(G)
			trQJ = 0
			for (bb in 1:ncol(Q)) {
				trQJ = trQJ + (Q[bb,] %*% J[,bb])
			}
			mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
			sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + 2*t(thetaHat)%*%(B - muVec[k,] %*% t(iVec)) %*% A + tr(thetaHat %*% t(thetaHat) %*% Omega) + gSigVec[i] *(1 + eta + trQJ + t(G) %*% C %*% G + tr(C %*% P - 2*C %*% t(hTilde) %*% invR %*% K))
			muLs[k, i] <- mu_L
			sigma2Ls[k, i] <- sigma2_L
			if(sigma2_L < 0) { 
				if(sigma2_L < -(10^-7)) {
					if(sigma2_L < -(10^-7)) {
						print(paste0("Negative variance at index ", k))
						print(paste0("Value of ", sigma2_L))
					}
				}
				sigma2Ls[k,i] = abs(sigma2_L)
			}
		}
	}
	samps = linkedSampling(fLocsPred, gLocsPred, muLs, sigma2Ls, fGamma, gGamma, numTest, gPred$mean)
	return(samps)
}
###
linkedConstantMultiple <- function(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta) {
	muLs <- matrix(c(0), nrow = numTest, ncol = g)
	sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
	sampsA <- array(0, dim = c(numTest, numTest, g))
	sampsSamples <- array(0, dim = c(numTest, 100, g))
	HzT <- matrix(1,m,1)
	hTilde <- matrix(1,m,1)
	hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% gOutput
	
	#calculate all progress locations
	progLocs <- ceil((1:10*(0.1))*(g*numTest))
	
	for (i in 1:g) {
		yT <- gOutput[,i]
		thetaHat <- rep(0, times = ncol(wT))
		betaHat <- hats[,i]

		sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde) #saving computation time since this doesn't change each iteration
		B <- matrix(0,nrow = d, ncol = m)
		Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
		C <- solve(t(hTilde) %*% invR %*% hTilde)
		A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)

		for(k in 1:numTest){
			currLoc <- ((i-1)*numTest) + k
			if(sum(currLoc == progLocs) == 1) { #print statement if 10/20/.../100% through calculations
				print(paste0(floor((currLoc/(g*numTest))*100), "% done"))
			}
			Sys.sleep(0.001) #can remove this when package is ready to go live probably
			flush.console() #can remove this when package is ready to go live probably
			hZ <- 1

			iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			K <- cbind(t(B), iVec%*%t(hZ))
			J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			tempOmega <- sigVec[k,]^2
			Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))

			P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
			G <- cbind(t(muVec[k,]), t(hZ))
			G <- t(G)
			trQJ = 0
			for (bb in 1:ncol(Q)) {
				trQJ = trQJ + (Q[bb,] %*% J[,bb])
			}
			mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
			sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + gSigVec[i]*(1 + eta + sigInvMat + trQJ - 2*tr(C * t(hTilde) %*% invR %*% iVec))
			muLs[k, i] <- mu_L
			sigma2Ls[k, i] <- sigma2_L
			if(sigma2_L < 0) { #add some sort of flag here
				print(paste0("Negative variance at index ", k))
				print(paste0("Value of ", sigma2_L))
				sigma2Ls[k,i] = abs(sigma2_L)
			}
		}
		tempSamps = linkedSampling(fLocsPred, gLocsPred, muLs[,i], sigma2Ls[,i], fGamma, gGamma, numTest, gPred$mean[,i])
		sampsA[,,i] = tempSamps$A
		sampsSamples[,,i] = tempSamps$samples
	}
	samps = list(A = sampsA, samples = sampsSamples)
	return(samps)
}
###
linkedLinearMultiple <- function(m, d, p, g, numTest, wT, zT, fLocsPred, fGamma, gLocsPred, gPred, gOutput, gGamma, muVec, sigVec, gSigVec, invR, eta) {
	muLs <- matrix(c(0), nrow = numTest, ncol = g)
	sigma2Ls <- matrix(c(0), nrow = numTest, ncol = g)
	sampsA <- array(0, dim = c(numTest, numTest, g))
	sampsSamples <- array(0, dim = c(numTest, 100, g))
	HzT <- cbind(zT,matrix(1,m,1))
	hTilde <- cbind(wT, HzT)
	hats <- solve(t(hTilde) %*% invR %*% hTilde)%*%t(hTilde) %*% invR %*% gOutput
	
	#calculate all progress locations
	progLocs <- ceil((1:10*(0.1))*(g*numTest))
	
	for (i in 1:g) {
		yT <- gOutput[,i]

		thetaHat <- hats[1:d,i]
		betaHat <- hats[(d+1):(d+p+1),i]

		sigInvMat <- solve(t(hTilde)%*%invR%*%hTilde) #saving computation time since this doesn't change each iteration
		B <- matrix(0,nrow = d, ncol = m)
		Q <- invR %*% hTilde %*% solve(t(hTilde) %*% invR %*% hTilde) %*% t(hTilde) %*% invR - invR
		C <- solve(t(hTilde) %*% invR %*% hTilde)

		A <- invR%*%(yT - wT %*% thetaHat - HzT %*% betaHat)
		for(k in 1:numTest){
			currLoc <- ((i-1)*numTest) + k
			if(sum(currLoc == progLocs) == 1) { #print statement if 10/20/.../100% through calculations
				print(paste0(floor((currLoc/(g*numTest))*100), "% done"))
			}
			Sys.sleep(0.001)
			flush.console()
			hZ <- c(gLocsPred[k,(d+1):(d+p)],1)

			iVec <- iFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			B <- bFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			K <- cbind(t(B), iVec%*%t(hZ))
			J <- jFun(gLocsPred[k,], as.matrix(wT), as.matrix(zT), gGamma, muVec[k,], sigVec[k,])
			tempOmega <- sigVec[k,]^2
			Omega <- matrix(diag(tempOmega), ncol = length(tempOmega))
			P <- blkdiag(Omega, matrix(0,nrow = length(hZ), ncol = length(hZ)))
			G <- cbind(t(muVec[k,]), t(hZ))
			G <- t(G)
			trQJ = 0
			for (bb in 1:ncol(Q)) {
				trQJ = trQJ + (Q[bb,] %*% J[,bb])
			}
			mu_L <- (t(muVec[k,]) %*% thetaHat) + (t(hZ) %*% betaHat) + (t(iVec) %*% A)
			sigma2_L <- t(A) %*% (J - iVec%*%t(iVec)) %*% A + 2*t(thetaHat)%*%(B - muVec[k,] %*% t(iVec)) %*% A + tr(thetaHat %*% t(thetaHat) %*% Omega) + gSigVec[i] *(1 + eta + trQJ + t(G) %*% C %*% G + tr(C %*% P - 2*C %*% t(hTilde) %*% invR %*% K))
			muLs[k, i] <- mu_L
			sigma2Ls[k, i] <- sigma2_L
			if(sigma2_L < 0) { 
				print(paste0("Negative variance at index ", k))
				print(paste0("Value of ", sigma2_L))
				sigma2Ls[k,i] = abs(sigma2_L)
			}
		}
		tempSamps = linkedSampling(fLocsPred, gLocsPred, muLs[,i], sigma2Ls[,i], fGamma, gGamma, numTest, gPred$mean[,i])
		sampsA[,,i] = tempSamps$A
		sampsSamples[,,i] = tempSamps$samples
	}
	samps = list(A = sampsA, samples = sampsSamples)
	return(samps)
}