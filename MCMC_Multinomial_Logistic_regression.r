	list(F=2, G= 2, H=2, I=2, J=2, K=5,
	X = structure(.Data = c(44, 22, 27, 56, 27, 51, 69, 57, 19, 45,  6, 77, 10, 90, 75, 37, 29, 65, 80, 51, 45, 12,  6, 27, 97, 10, 33, 46, 77, 22, 42,  9, 98, 98, 23, 59, 78, 60, 14, 41, 55, 26, 30, 42, 23, 19, 62, 13, 36, 11, 57, 97, 45, 49, 68, 98, 29, 86, 26, 94, 56, 37, 37, 26, 18, 19, 30, 36, 32, 51, 70, 91, 46, 22, 40, 38,  3, 90, 77, 55, 11, 11,  3, 52, 16, 52, 99, 67, 73, 96, 25, 52, 12, 71, 71, 43, 92, 76, 76, 96, 23, 21, 49, 56, 51, 62, 67, 25, 86, 80, 86, 62, 27, 26,  3, 54, 63, 21, 48, 69,  2,  2, 27, 61, 33, 84,  3, 78, 87, 70, 67, 68, 14, 12, 97, 88, 47, 77, 79, 48, 71,  9, 30, 34, 93, 54, 18, 20, 57, 60, 73, 22, 77, 28, 34, 84, 43, 95, 54, 28), .Dim = c(2, 2, 2, 2, 2, 5)))


	list(alpha = c(NA, 0, 0, 0, 0),
		beta = structure(.Data = c(NA, NA, NA, NA, NA,
     										       NA, 0, 0, 0, 0), .Dim = c(2, 5)),
		gamma = structure(.Data = c(NA, NA, NA, NA, NA,
     										       NA, 0, 0, 0, 0), .Dim = c(2, 5)),
		delta = structure(.Data = c(NA, NA, NA, NA, NA,
     										       NA, 0, 0, 0, 0), .Dim = c(2, 5)),
		epsilon = structure(.Data = c(NA, NA, NA, NA, NA,
     										       NA, 0, 0, 0, 0), .Dim = c(2, 5)),
		zita = structure(.Data = c(NA, NA, NA, NA, NA,
     										       NA, 0, 0, 0, 0), .Dim = c(2, 5)),
		lambda = structure(.Data = c(0, 0,
									 0, 0,
									 0, 0,
									 0, 0,
						
	model
	{

	#  PRIORS
		alpha[1] <- 0;       # zero contrast for baseline food
		for (k in 2 : K) { 
			alpha[k] ~ dnorm(0, 0.00001) # vague priors
		} 


	# Loop around first security measure(1rst):                                        1 => beta > f dans F
		for (k in 1 : K){  
			beta[1, k] <- 0 
		} # corner-point contrast with first lake 
		for (f in 2 : F) {     
			beta[f, 1] <- 0 ;  # zero contrast for baseline food
			for (k in 2 : K){  
				beta[f, k] ~ dnorm(0, 0.00001) # vague priors
			} 
		}
	# Loop around second security measure(2nd):                                       2 => gamma > g dans G
		for (k in 1 : K){  
			gamma[1, k] <- 0 
		} # corner-point contrast with first lake 
		for (g in 2 : G) {     
			gamma[g, 1] <- 0 ;  # zero contrast for baseline food
			for (k in 2 : K){  
				gamma[g, k] ~ dnorm(0, 0.00001) # vague priors
			} 
		}
	# Loop around third security mesure (3rd) :                                        3 => delta > h dans H
		for (k in 1 : K){  
			delta[1, k] <- 0 # corner-point contrast with first size 
		}  
		for (h in 2 : H){     
			delta[h, 1] <- 0 ;  # zero contrast for baseline food
			for (k in 2 : K){ 
				delta[h, k] ~ dnorm(0, 0.00001) # vague priors
			} 
		}

	# Loop around fourth security mesure (4nd) :                                        4 => epsilon > i dans I
		for (k in 1 : K){  
			epsilon[1, k] <- 0 # corner-point contrast with first size 
		}  
		for (i in 2 : I){     
			epsilon[i, 1] <- 0 ;  # zero contrast for baseline food
			for (k in 2 : K){ 
				epsilon[i, k] ~ dnorm(0, 0.00001) # vague priors
			} 
		}

	# Loop around fithdt security mesure (5nd) :                                        5 => zita > j dans J
		for (k in 1 : K){  
			zita[1, k] <- 0 # corner-point contrast with first size 
		}  
		for (j in 2 : J) {     
			zita[j, 1] <- 0 ;  # zero contrast for baseline food
			for ( k in 2 : K){ 
				zita[j, k] ~ dnorm(0, 0.00001) # vague priors
			} 
		}

	# LIKELIHOOD	
		for (f in 1 : F) {     # loop around lakes
			for (g in 1 : G) {     # loop around sizes
				for (h in 1 : H) {
					for (i in 1 : I) {
						for (j in 1 : J) {
							
	# Multinomial response
	      X[f, g, h, i, j ,1 : K] ~ dmulti( p[f, g, h, i, j, 1 : K] , n[f, g, h, i, j]  )
	      n[f, g, h, i, j] <- sum(X[f, g, h, i, j, ])
	      for (k in 1 : K) {     # loop around foods
	         p[f, g, h, i, j, k]        <- phi[f, g, h, i, j, k] / sum(phi[f, g, h, i, j, ])
	         log(phi[f, g, h, i ,j, k]) <- alpha[k] + beta[f, k]  + gamma[g, k] + delta[h, k]  + epsilon[i, k]  + zita[j, k] 
	        }
					}  
				}
			}  
		}
	} 
	# TRANSFORM OUTPUT TO ENABLE COMPARISON 
	#  WITH AGRESTI'S RESULTS

		for (k in 1 : K) {     # loop around foods
			for (f in 1 : F) {     # loop around lakes
				msre_1[f, k] <- beta[f, k] - mean(beta[, k]);   # sum to zero constraint
			}
			for (g in 1 : G) {     # loop around sizes
				msre_2[g, k] <- gamma[g, k] - mean(gamma[, k]); # sum to zero constraint
			}
			for (h in 1 : H) {     # loop around sizes
				msre_3[h, k] <- delta[h, k] - mean(delta[, k]); # sum to zero constraint
			}
			for (i in 1 : I) {     # loop around sizes
				msre_4[i, k] <- epsilon[i, k] - mean(epsilon[, k]); # sum to zero constraint
			}
			for (j in 1 : J) {     # loop around sizes
				msre_5[j, k] <- zita[j, k] - mean(zita[, k]); # sum to zero constraint
			}
		}
	} 