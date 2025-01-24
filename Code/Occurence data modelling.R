#-------------------------------------------------------------------------
# program reads data on reported distributional characteristics of intake
# concentrations on T2/HT2 in food
# examples of charactertistics are overall mean, mean of measured values,
# maximum value, etc.
# several fit-functions (cases) are defined that are used to fit a log-normal distribution
# on the characteristics; the functions differ in the characteristics (information) being
# available for a study (row in database)
#-------------------------------------------------------------------------

rm(list = ls())
library(openxlsx)
library(nlme)
eps			<- 1e-3					# small constant
minlogmu		<- -10
LODQimputind	<- 2
f.LODQimput <- function(LOD2, LOQ2) {
	if (LODQimputind == 1) { return(exp(LOD2)) }
	if (LODQimputind == 2) { return(max(exp(LOD2), .5 * exp(LOQ2))) }
	}
currdate		<- paste("06122024", LODQimputind, sep = "")

#-------------------------------------------------------------------------
# calculation of boundary values
#-------------------------------------------------------------------------

f.makeUB1 <- function(hc) { return(as.numeric(strsplit(hc, "<")[[1]][2])) }
f.makeUB2 <- function(i, hdat5) {
	hc1	<- hdat5[i]
	if (is.na(hc1)) { hres <- NA
	} else {
		if (nchar(hc1) > 1) {
			while (substr(hc1, 1, 1) == " ") { hc1 <- substr(hc1, 2, nchar(hc1)) }}
		if (hc1 == "<LOQ") { hres <- LOQ[i] } else { hres <- NA }
	}
	return(hres) }

#-------------------------------------------------------------------------
# calculation of geometric mean values
#-------------------------------------------------------------------------

f.calcgeomean <- function(ch1, hsplit) {
	ch2	<-strsplit(ch1,hsplit)[[1]]
	x 	<- ch1
	if (length(ch2) > 1) {
		hb<-prod(!is.na(as.numeric(ch2)))
		if (hb) { x <- exp(mean(log(as.numeric(ch2)))) } else { x <- NA }
	}
	return(x) }

#-------------------------------------------------------------------------
# reading excel-file with info on distributional charactertistics
# rows: experiments, columns: distributional charactertistics
#-------------------------------------------------------------------------

hfile			<- "Occurrence data.xlsx"
hdat			<- read.xlsx(hfile, 1)
hnrow			<- dim(hdat)[1]

#-------------------------------------------------------------------------
# impute geometric mean values
#-------------------------------------------------------------------------

for (i1 in 1:hnrow) { for (i2 in 5:11) { 
	hdat[i1, i2] <- f.calcgeomean(hdat[i1, i2], "-")
	hdat[i1, i2] <- f.calcgeomean(hdat[i1, i2], " or ")
}}

hdat[117, 5]	<- hdat[121,5] <- NA

#-------------------------------------------------------------------------
# identify columns
#-------------------------------------------------------------------------

hnrow			<- dim(hdat)[1]				# # rows
subst			<- hdat[[1]]				# mycotoxin
rawfood		<- hdat[[2]]				# occurrence in raw food
procfood		<- hdat[[3]]				# occurrence in processed food
medlevel		<- as.numeric(hdat[[4]])		# median level
meanlevel		<- as.numeric(hdat[[5]])		# average level
meanlevelbelowLOQ	<- hdat[[5]] == "<LOQ"			# average level below LOQ
meanlevelbelowLOQ[is.na(meanlevelbelowLOQ)] <- FALSE
meanlevelpos	<- as.numeric(hdat[[6]])		# average level (in positive samples only)
maxlevel		<- as.numeric(hdat[[7]])		# max level
poscase		<- as.numeric(hdat[[8]]) / 100	# positive fraction of sample
LOQLODind		<- hdat[[9]]				# <LOD/LOQ
LOD			<- as.numeric(hdat[[10]])		# LOD level
LOQ			<- as.numeric(hdat[[11]])		# LOQ level
meanlevelUB1	<- unlist(lapply(as.list(hdat[[5]]), f.makeUB1)) # UB of average level
meanlevelUB2	<- unlist(lapply(as.list(1:hnrow), f.makeUB2, hdat5 = hdat[[5]]))
meanlevelUB		<- ifelse(!is.na(meanlevelUB1), meanlevelUB1, meanlevelUB2)
method		<- hdat[[12]]				# method of analysis
country		<- hdat[[13]]				# country
continent		<- hdat[[14]]				# contintent
refer			<- hdat[[16]]				# reference

meanlevelposUB	<- unlist(lapply(as.list(1:hnrow), f.makeUB2, hdat5 = hdat[[6]]))
meanlevelposUBind	<- which(!is.na(meanlevelposUB))
meanlevelpos[meanlevelposUBind] <- as.numeric(meanlevelposUB[meanlevelposUBind])

poscase[(LOQLODind == "<LOD") & is.na(poscase)] <- 0

#-------------------------------------------------------------------------
# reading excel-file with info on sampe sizes of these experiments
#-------------------------------------------------------------------------

hfile1		<- "Overview occurrence T2HT2.xlsx"
hdat1			<- read.xlsx(hfile1, 1)
hnrow1		<- dim(hdat1)[1]				# # rows
refer1		<- hdat1[[17]]				# reference
nsample1		<- as.numeric(hdat1[[4]])		# # samples

hest			<- array(0, dim = c(hnrow, 9))	# output array
hest[ , 1]		<- -1
hcorrmusig		<- rep(NA, hnrow)

#-------------------------------------------------------------------------
# Definitions of all fit-functions
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# fit function for case 0: everything known
#-------------------------------------------------------------------------

f.fun0 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax <- pnorm(maxlevel0, mu, sig) 													# cumulative probability of being under maximum
	fmax <- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 <- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 <- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +						# calculated mean value
		FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)		
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +					# calculated variance
		FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8)

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * exp(LOD0)) / (FLODcomp)	# calculated mean value over >LOD
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * exp(LOD0)^2) / (FLODcomp) - hmeanpos^2, 1E-8)	# calculated variance

	Fmed 	<- pnorm(medlevel0, mu, sig)													# probability of median level

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample) 					# term related to mean
	term4 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample)) 			# term related to mean
	term5 <- -.5 * nsample * log(Fmed) - .5 * nsample * log(1 - Fmed) 							# term related to median

	return(term1 + term2 + term3 + term4 + term5)
}

#-------------------------------------------------------------------------
# fit function for case 1: known LOD&LOQ, max value, # pos values, and overall mean value
#-------------------------------------------------------------------------

f.fun1 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax <- pnorm(maxlevel0, mu, sig) 													# cumulative probability of being under maximum
	fmax <- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 <- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 <- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +						# calculated mean value
		FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)		
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +					# calculated variance
		FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8)

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample) 						# term related to mean
	
	return(term1 + term2 + term3)
}

#-------------------------------------------------------------------------
# fit function for case 12: everything except LOD/LOQ
#-------------------------------------------------------------------------

f.fun12 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	LOD0 	<- par1[3]																# LOD-value
	
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax <- pnorm(maxlevel0, mu, sig) 													# cumulative probability of being under maximum
	fmax <- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 <- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 <- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +						# calculated mean value
		FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)	
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +					# calculated variance
		FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8)

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)) / (FLODcomp)	# calculated mean value over >LOD
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)^2) / (FLODcomp) - hmeanpos^2, 1E-8)	# calculated variance
						
	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample) 						# term related to mean
	term4 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample)) 			# term related to mean

	return(term1 + term2 + term3 + term4)
}

#-------------------------------------------------------------------------
# fit function for case 11: known LOD&LOQ, max value, # pos values
#-------------------------------------------------------------------------

f.fun11 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax <- pnorm(maxlevel0, mu, sig) 													# cumulative probability of being under maximum
	fmax <- dnorm(maxlevel0, mu, sig) 													# related probability density
	
	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	
	return(term1 + term2)
}

#-------------------------------------------------------------------------
# fit function for case 2: unknown LOQ&LOD, max value, # pos values, and overall mean value
#-------------------------------------------------------------------------

f.fun2 <- function(par1) {
	mu 	<- par1[1]																
	sig 	<- exp(par1[2])
	LOD0 	<- par1[3]																# LOD-value
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	LOQ0	<- log(2 * exp(LOD0))
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 	<- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 	<- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +						# calculated mean value
		FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +					# calculated variance
		FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8)

	term1 <<- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term2 <<- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <<- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample)		 				# term related to mean
	
	return(term1 + term2 + term3)
}

#-------------------------------------------------------------------------
# fit function for case 7: known LOQ, max value, # pos values, and UB overall mean value; not used
#-------------------------------------------------------------------------

f.fun7 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])

	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density

	p1 	<- pnorm(meanlevelUB0, mu, sig / sqrt(nsample))

	f.x1 	<- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 	<- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD  * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)						# calculated mean value
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD  * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8) 		# calculated variance

	term1 <- -nsample * (1 - poscase0) * log(FLOQ) - poscase0 * nsample * log(FLOQcomp) 				# term related to LOD
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- -log(p1)																# term related to reported upper boundary of mean value

	return(term1 + term2 + term3)
}

#-------------------------------------------------------------------------
# fit function for case 3: known LOD&LOQ, max value, # pos values, and mean value over >LOD values, overall mean value missing
#-------------------------------------------------------------------------

f.fun3 <-function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 	<- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 	<- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)) / (FLODcomp)	# calculated mean value over >LOD
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)^2) / (FLODcomp) - hmeanpos^2, 1E-8)	# calculated variance
																
	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOQ
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample)) 			# term related to mean
	
	return(term1 + term2 + term3)
}


#-------------------------------------------------------------------------
# fit function for case 4: known LOD$LOQ, max value, # pos values, and mean value over >LOD values, overall mean value <LOQ
#-------------------------------------------------------------------------

f.fun4 <-function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	muexp <- exp(mu + .5 * sig^2)
	sigexp <- muexp * sqrt(exp(sig^2) - 1)

	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 	<- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 	<- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)) / (FLOQcomp)	# calculated mean value over >LOD
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)^2)/ (FLOQcomp) - hmeanpos^2, 1E-8) # calculated variance

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOQ
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample)) 			# term related to mean over >LOQ
	term4 <- pnorm(exp(LOD0), muexp, sigexp / sqrt(nsample))									# term related to overall mean

	return(term1 + term2 + term3 + term4)
}

#-------------------------------------------------------------------------
# fit function for case 9: unknown LOD&LOQ, max value, # pos values, and mean value over >LOD values
#-------------------------------------------------------------------------

f.fun9 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])

	LOQD 	<- par1[3]
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	LOQ0	<- log(2 * exp(LOD0))
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 	<- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 	<- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)) / (FLODcomp)	# calculated mean value
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)^2) / (FLODcomp) - hmeanpos^2, 1E-8) # calculated variance

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOQ
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample)) 			# term related to mean

	return(term1 + term2 + term3)
}

#-------------------------------------------------------------------------
# fit function for case 5: known LOD&LOQ, max value, # pos values, mean value over all and over >LOD values
#-------------------------------------------------------------------------

f.fun5 <-function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])

	FLOD <<- max(pnorm(LOD0, mu, sig), 1E-8) 												# probability of <LOD
	FLODcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related density

	f.x1 	<- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 	<- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)						# calculated overall mean value
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8) 		# calculated overall variance

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)) / (FLODcomp)	# calculated mean value over positive values
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)^2) / (FLODcomp) - hmeanpos^2, 1E-8) # calculated variance over positive values
																
	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOQ
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample)) 			# term related to mean
	term4 <- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample) 					# term related to mean

	return(term1 + term2 + term3 + term4)
}

#-------------------------------------------------------------------------
# fit function for case 6: unknown LOD&LOQ, max value, # pos values, mean value over all and over >LOQ values
#-------------------------------------------------------------------------

f.fun6 <-function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	LOD0 	<- par1[3]
	LOQ0	<- log(2 * exp(LOD0))

	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of not being positive
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of being positive
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 	<- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 	<- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)						# calculated mean value
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8) 		# calculated variance

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)) / (FLODcomp) # calculated mean value over positive values
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)^2) / (FLODcomp) - hmeanpos^2, 1E-8) # calculated variance over positive values

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOQ
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample)) 			# term related to overall mean
	term4 <- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample) 					# term related to mean over positive values

	return(term1 + term2 + term3 + term4)
}

#-------------------------------------------------------------------------
# fit function for case 8: known LOD&LOQ, max value, # pos values, and overall median value
#-------------------------------------------------------------------------

f.fun8 <-function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])

	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density
	Fmed 	<- pnorm(medlevel0, mu, sig)													# probability of median level

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOQ
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- -.5 * nsample * log(Fmed) - .5 * nsample * log(1 - Fmed) 							# term related to median

	return(term1 + term2 + term3)
}

#-------------------------------------------------------------------------
# case 10: only positive rate
#-------------------------------------------------------------------------

f.fun10 <- function(par1) {
	mu	<- par1[1]

	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of not being positive
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of being positive
	
	Fmax 	<- pnorm(maxlevel0, mu, sig) 													# probability of being under maximum
	fmax 	<- dnorm(maxlevel0, mu, sig) 													# related probability density

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	
	return(term1) }

#-------------------------------------------------------------------------
# case 13: known max-level, fraction positives, LOQ, and mean level
#-------------------------------------------------------------------------

f.fun13 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	Fmax <- pnorm(maxlevel0, mu, sig) 													# cumulative probability of being under maximum
	fmax <- dnorm(maxlevel0, mu, sig) 													# related probability density

	f.x1 <- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 <- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)						# calculated mean value
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8) 		# calculated variance

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term2 <- -(nsample - 1) * log(Fmax) - log(fmax) 										# term related to maximum
	term3 <- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample) 					# term related to mean
	
	return(term1 + term2 + term3)
}

#-------------------------------------------------------------------------
# case 14: known fraction positives, LOQ, and mean level
#-------------------------------------------------------------------------

f.fun14 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	f.x1 <- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 <- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanall <<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * .5 * exp(LOD0) + FLODQ * f.LODQimput(LOD0, LOQ0)						# calculated mean value
	hvarall <<- max(integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val +
				FLOD * (.5 * exp(LOD0))^2 + FLODQ * f.LODQimput(LOD0, LOQ0)^2 - hmeanall^2, 1E-8) 		# calculated variance

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term3 <- .5 * log(hvarall) + .5 * (mean0all - hmeanall)^2 / (hvarall / nsample) 					# term related to mean
	
	return(term1 + term3)
}

#-------------------------------------------------------------------------
# case 15: known fraction positive, mean value over positives, and LOQ
#-------------------------------------------------------------------------

f.fun15 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- exp(par1[2])
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ

	f.x1 <- function(x, mu1, sig1) { exp(x) * dnorm(x, mu1, sig1) } 								# function to calculate mean value
	f.x2 <- function(x, mu1, sig1) { exp(x)^2 * dnorm(x, mu1, sig1) } 							# function to calculate mean quadratic value

	hmeanpos <<- (integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)) / (FLODcomp) # calculated mean value over positive values
	hvarpos <<- max((integrate(f.x2, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val + FLODQ * f.LODQimput(LOD0, LOQ0)^2) / (FLODcomp) - hmeanpos^2, 1E-8) # calculated variance over positive values

	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	term3 <- .5 * log(hvarpos) + .5 * (mean0pos - hmeanpos)^2 / (hvarpos / (poscase0 * nsample))			# term related to mean
	
	return(term1 + term3)
}
#-------------------------------------------------------------------------
# case 17: known LOQ, and some positive cases
#-------------------------------------------------------------------------

f.fun17 <- function(par1) {
	mu 	<- par1[1]
	sig 	<- sigopt
	
	FLOD 	<<- max(pnorm(LOD0, mu, sig), 1E-8) 											# probability of <LOD
	FLODcomp <- max(1 - FLOD, 1E-8)													# probability of >LOD
	FLOQ 	<<- max(pnorm(LOQ0, mu, sig), 1E-8) 											# probability of <LOQ
	FLOQcomp <- max(1 - FLOQ, 1E-8)													# probability of >LOQ
	FLODQ <- FLOQ - FLOD															# probability of LOD<.<LOQ	
					
	term1 <- -nsample * (1 - poscase0) * log(FLOD) - poscase0 * nsample * log(FLODcomp) 				# term related to LOD
	
	return(term1)
}

#-------------------------------------------------------------------------
# case 16: known LOQ, and no positive cases, e.g. Triticale
#-------------------------------------------------------------------------

f.fun16 <- function(hsig) {
	pest	<- pnorm(LOD0, mu3, hsig)
	return(sum((pest - pest0)^2)) }

f.calcvar <- function(f.fun, hpar, eps) {
	f0	<- f.fun(hpar)
	f1	<- f.fun((1 + eps) * hpar)
	f2	<- f.fun((1 - eps) * hpar)
	hvar	<- (f1 - 2 * f0 + f2) / (hpar * eps)^2
	return(1/hvar) }

#-------------------------------------------------------------------------
# Imputations of missing LOD or LOQ values
#-------------------------------------------------------------------------

hind		<- !is.na(LOD) & is.na(LOQ)
LOQ[hind]	<- LOD[hind]
hind		<- is.na(LOD) & !is.na(LOQ)
LOD[hind]	<- LOQ[hind]

#-------------------------------------------------------------------------
# Definitions of cases, selecting rows that meet the case criteria, 
# applying the fit-function, some post-calculations, and stoarge of results in result-file
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# case 0: available: overall mean level, max level, # pos cases, LOQ
#-------------------------------------------------------------------------

hind0	<- which(!is.na(medlevel) & !is.na(meanlevel) & !is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind0) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	mean0all 	<- meanlevel[i]						# untransformed mean value
	logmean0 	<- log(mean0all)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4)	# starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])					# log-transformed max value
	mean0pos 	<- meanlevelpos[i]					# log-transformed mean value
	medlevel0 	<- log(medlevel[i])				# log-transformed median value

	hres	<- optim(par0, f.fun0, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov		<- solve(hres$hessian)				# co-variance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2]	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		hres1 	<- c(0, hpar, hse, 1 - FLOD, hmeanall)		
		print(hres1)
		hest[i, 1:7] <- hres1 						# estimation results stored in output array
	}
}}



#-------------------------------------------------------------------------
# case 1: available: overall mean level, max level, # pos cases, LOQ
#-------------------------------------------------------------------------

hind1	<- which(is.na(medlevel) & !is.na(meanlevel) & is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind1) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	mean0all 	<- meanlevel[i]						# untransformed mean value
	logmean0 	<- log(mean0all)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4)	# starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])					# log-transformed max value
		
	hres	<- optim(par0, f.fun1, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov		<- solve(hres$hessian)				# co-variance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2]	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		hres1 	<- c(1, hpar, hse, 1 - FLOD, hmeanall)		
		print(hres1)
		hest[i, 1:7] <- hres1 						# estimation results stored in output array
	}
}}

#-------------------------------------------------------------------------
# case 2: available: overall mean level, max level, # pos cases, unknown LOQ
#-------------------------------------------------------------------------

hind2	<- which(is.na(medlevel) & !is.na(meanlevel) & is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & is.na(LOQ) & !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind2) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	mean0all 	<- meanlevel[i]						# log-transformed mean value
	logmean0 	<- log(mean0all)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4, logmean0 / 2) # starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])					# log-transformed max value
	
	hres	<- optim(par0, f.fun2, hessian = TRUE, control = list(maxit = 1000)) # optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {  						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian)				# parameter estimation covariance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2] 	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		hres1 	<- c(2, hpar[1:2], hse[1:2], 1 - FLOD, hmeanall)		
		print(hres1)
		hest[i, 1:7] <- hres1 						# estimation results stored in output array
	}
}}

#-------------------------------------------------------------------------
# case 3: available: mean level > LOQ, max level, # pos cases, LOQ
#-------------------------------------------------------------------------

hind3	<- which(is.na(medlevel) & is.na(meanlevel) & !is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & !meanlevelbelowLOQ & (hest[ , 1] == -1))

if (TRUE) { for (i in hind3) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	LOD0 		<- log(LOD[i])						# log-transformed LOD-value 
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value 
	mean0pos 	<- meanlevelpos[i]					# log-transformed mean value
	logmean0 	<- log(mean0pos)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4) # starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])					# log-transformed max value

	hres 	<- optim(par0, f.fun3, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian)				# parameter estimation covariance matrix
		if (min(diag(hcov)) > 0) {					# acceptable covariance matrix
			hcov 		<- solve(hres$hessian)			# parameter estimation covariance matrix
			hse 		<- sqrt(diag(hcov ))			# calculated standard errors of parameter estimates
			hcorrmusig[i] <- hcov[1, 2] / prod(hse)		# correlation between estimated mu and sig
			hse[2] 	<- hpar[2] * hse[2]			# adjustment for sigma-parameter
			hest[i, c(1:6, 8)] <- c(3, hpar, hse, 1 - FLOD, hmeanpos)			
		}
	}
}}

#-------------------------------------------------------------------------
# case 4: available: mean level > LOQ, max level, # pos cases, LOQ, overall mean<LOQ
#-------------------------------------------------------------------------

hind4 <- which(is.na(medlevel) & is.na(meanlevel) & !is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & meanlevelbelowLOQ & (hest[ , 1] == -1))

if (TRUE) { for (i in hind4) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	LOD0 		<- log(LOD[i])					# log-transformed LOQ-value 
	LOQ0 		<- log(LOQ[i])					# log-transformed LOQ-value 
	mean0pos 	<- meanlevelpos[i]				# log-transformed mean value
	logmean0 	<- log(mean0pos)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4) # starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])				# log-transformed max value

	hres 	<- optim(par0, f.fun4, hessian = TRUE) 		# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {					# acceptable model fit
		hpar[2] 	<- exp(hpar[2])				# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian)			# parameter estimation covariance matrix
		if (min(diag(hcov)) > 0) {				# acceptable covariance matrix
			hcov 		<- solve(hres$hessian)		# parameter estimation covariance matrix
			hse 		<- sqrt(diag(hcov ))		# calculated standard errors of parameter estimates
			hcorrmusig[i] <- hcov[1, 2] / prod(hse)		# correlation between estimated mu and sig
			hse[2] 	<- hpar[2] * hse[2]		# adjustment for sigma-parameter
			hest[i, c(1:6, 8)] <- c(4, hpar, hse, 1 - FLOD, hmeanpos) 			
		}
	}
}}

#-------------------------------------------------------------------------
# case 5: available: mean level overall and > LOQ, max level, # pos cases, LOQ
#-------------------------------------------------------------------------

hind5	<- which(is.na(medlevel) & !is.na(meanlevel) & !is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & !meanlevelbelowLOQ & (hest[ , 1] == -1))

if (TRUE) { for (i in hind5) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	LOD0 		<- log(LOD[i])					# log-transformed LOD-value		
	LOQ0 		<- log(LOQ[i])					# log-transformed LOQ-value		
	mean0pos 	<- meanlevelpos[i]				# log-transformed mean value over positive values
	mean0all 	<- meanlevel[i]					# log-transformed overall mean value
	logmean0 	<- log(mean0all)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4) # starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])				# log-transformed max value

	hres 	<- optim(par0, f.fun5, hessian = TRUE) 		# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {					# acceptable model fit
		hpar[2] 	<- exp(hpar[2])				# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian)			# parameter estimation covariance matrix
		hse		<- sqrt(diag(hcov))			# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)		# correlation between estimated mu and sig
		hse[2] 	<- hpar[2] * hse[2]			# adjustment for sigma-parameter
		hest[i, 1:8] <- c(5, hpar, hse, 1 - FLOD, hmeanall, hmeanpos) 		
	}
}}

#-------------------------------------------------------------------------
# case 6: available: mean level overall and > LOQ, max level, # pos cases, unknown LOQ
#-------------------------------------------------------------------------

hind6	<- which(is.na(medlevel) & !is.na(meanlevel) & !is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & is.na(LOQ) & !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind6) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	mean0pos 	<- meanlevelpos[i]				# log-transformed mean value over positive values		
	mean0 	<- meanlevel[i]					# log-transformed overall mean value
	logmean0 	<- log(mean0all)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4, logmean0 / 2) # starting vector for optimization
	maxlevel0 <- log(maxlevel[i])					# log-transformed max value	

	hres 	<- optim(par0, f.fun6, hessian = TRUE) 		# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {					# acceptable model fit					
		hpar[2] 	<- exp(hpar[2])				# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian)			# parameter estimation covariance matrix
		hse		<- sqrt(diag(hcov))			# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)		# correlation between estimated mu and sig
		hse[2] 	<- hpar[2] * hse[2]			# adjustment for sigma-parameter
		hest[i, 1:8] <- c(6, hpar[1:2], hse[1:2], 1 - FLOD, hmeanall, hmeanpos) 		
	}
}}

#-------------------------------------------------------------------------
# case 7: available: UB mean level overall, max level, # pos cases, LOQ
#-------------------------------------------------------------------------

hind7	<- which(is.na(medlevel) & is.na(meanlevel) & is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & !is.na(meanlevelUB) & (hest[ , 1] == -1))


if (TRUE) { for (i in hind7) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	LOD0 		<- log(LOD[i])					# log-transformed LOD-value	
	LOQ0 		<- log(LOQ[i])					# log-transformed LOQ-value	
	meanlevelUB0 <- log(meanlevelUB[i])				# log-transformed UB mean level 
	par0 		<- c(meanlevelUB0 + eps, abs(meanlevelUB0 + eps) / 4) # starting vector for optimization
	maxlevel0	 <- log(maxlevel[i])				# log-transformed max value

	hres 	<- optim(par0, f.fun7, hessian = TRUE) 		# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])				# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian)			# parameter estimation covariance matrix
		hse		<- sqrt(diag(hcov ))			# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)		# correlation between estimated mu and sig
		hse[2] 	<- hpar[2] * hse[2]			# adjustment for sigma-parameter
		hest[i, 1:7] <- c(7, hpar, hse, 1 - FLOD, hmeanall) 
		
	}
}}

#-------------------------------------------------------------------------
# case 8: available: median level, max level, # pos cases, LOQ
#-------------------------------------------------------------------------

hind8	<- which(!is.na(medlevel) & is.na(meanlevel) & is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind8) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	LOD0 		<- log(LOD[i])					# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])					# log-transformed LOQ-value
	medlevel0 	<- log(medlevel[i])				# log-transformed median value
	par0 		<- c(medlevel0 + eps, abs(medlevel0 + eps) / 4) # starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])				# log-transformed max value	

	hres <- optim(par0, f.fun8, hessian = TRUE) 		# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {					# acceptable model fit
		hpar[2] 	<- exp(hpar[2])				# back-transformation of sigma-parameter
		mu 		<- hpar[1]
		sig 		<- hpar[2]
		f.x1 <- function(x, mu1, sig1) { x * dnorm(x, mu1, sig1) }	# function to calculate mean value
		hmeanall 	<- integrate(f.x1, LOQ0, mu + 6 * sig, mu1 = mu, sig1 = sig)$val 
		hmeanpos 	<- hmeanall / (1 - pnorm(LOQ0, mu, sig))		# calculated mean over positieve values
		hcov		<- solve(hres$hessian)			# co-variance matrix
		hse		<- sqrt(diag(hcov))			# standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)		# correlation between estimated mu and sig
		hse[2] 	<- hpar[2] * hse[2]			# adjustment for sigma-parameter
		hest[i, 1:9] <- c(8, hpar, hse, 1 - FLOD, hmeanall, hmeanpos, exp(hpar[1]))
	}
}}

#-------------------------------------------------------------------------
# case 9: available: mean level > LOD, max level, # pos cases,
#-------------------------------------------------------------------------

hind9 <- which(is.na(medlevel) & is.na(meanlevel) & !is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & is.na(LOQ) & !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind9) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	mean0pos 	<- meanlevelpos[i]				# log-transformed mean value
	logmean0 	<- log(mean0pos)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4, logmean0 / 2) # starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])				# log-transformed max value

	hres <- optim(par0, f.fun9, hessian = TRUE) 		# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {					# acceptable model fit
		hpar[2] 	<- exp(hpar[2])				# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian[1:2,1:2])	# parameter estimation covariance matrix
		if (min(diag(hcov)) > 0) {				# acceptable covariance matrix
			hse 		<- sqrt(diag(hcov))		# standard errors of parameter estimates
			hcorrmusig[i] <- hcov[1, 2] / prod(hse)	# correlation between estimated mu and sig
			hse[2] 	<- hpar[2] * hse[2]		# adjustment for sigma-parameter
			hest[i, c(1:6, 8)] <- c(9, hpar[1:2], hse, 1 - FLOD, hmeanpos) 
		}
	}
}}

#-------------------------------------------------------------------------
# case 10: only positive rate
#-------------------------------------------------------------------------

hind10 <- which(is.na(medlevel) & is.na(meanlevel) & is.na(meanlevelpos) & is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & (LOQLODind == "<LOD") & (hest[ , 1] == -1))

colfood 	<- 2
food0		<- intersect(unique(hdat[hest[ , 1] == 1, colfood]), unique(hdat[hind9, colfood]))
nfood0	<- length(food0)
if (nfood0 == 0) {
	colfood 	<- 3
	food0		<- intersect(unique(hdat[hest[ , 1] == 1, colfood]), unique(hdat[hind9, colfood]))
	nfood0	<- length(food0) }

for (i0 in seq(length = nfood0)) {
	food		<- food0[i0]
	sigall	<- hest[which((hest[ , 1] == 1) & (hdat[ , 2] == food)), 3]
	sig		<- mean(sigall)
	sigse		<- sqrt(var(sigall))
	hind9A	<- hind9[hdat[hind9, 2] == food]
	for (i in hind9A) {
		LOD0		<- log(LOD[i])
		nsample	<- nsample1[i]
		poscase0	<- poscase[i] 
		par0 		<- c(.5, 2) * LOD0 / 2
		hres	 	<- optimize(f.fun10, par0) # optimization
		paropt	<- hres$minimum
		hse		<- sqrt(f.calcvar(f.fun10, paropt, .1))
		hest[i, 1:5] <- c(10, paropt, sig, hse, sigse)
	}
}
hind10	<- which(hest[ , 1] == 10)

#-------------------------------------------------------------------------
# case 11: max, positive
#-------------------------------------------------------------------------

hind11 <- which(is.na(medlevel) & is.na(meanlevel) & is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ)& !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind11) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	mean0 	<- LOQ[i] / 2						# untransformed mean value
	logmean0 	<- log(mean0)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4)	# starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])					# log-transformed max value
		
	hres	<- optim(par0, f.fun11, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov		<- solve(hres$hessian)				# co-variance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2]	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		hres1 	<- c(11, hpar, hse, 1 - FLOD, hmeanall)		
		print(hres1)
		hest[i, 1:7] <- hres1 						# estimation results stored in output array
	}
}}

#-------------------------------------------------------------------------
# case 12: everything except (borders of) LOD/LOQ
#-------------------------------------------------------------------------

hind12 <- which(!is.na(medlevel) & !is.na(meanlevel) & !is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & is.na(LOQ)& !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind12) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]

	mean0pos 	<- meanlevelpos[i]					# log-transformed mean value over positive values		
	mean0all 	<- meanlevel[i]						# log-transformed overall mean value
	logmean0 	<- log(mean0all)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4, logmean0 / 2) # starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])					# log-transformed max value	

	hres 	<- optim(par0, f.fun12, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit					
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov 		<- solve(hres$hessian)				# parameter estimation covariance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2] 	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		hest[i, 1:8] <- c(12, hpar[1:2], hse[1:2], 1 - FLOD, hmeanall, hmeanpos) 		
	}
}}

#-------------------------------------------------------------------------
# case 13: known max-level, fraction positives, LOQ, and mean level
#-------------------------------------------------------------------------

hind13 <- which(is.na(medlevel) & is.na(meanlevel) & is.na(meanlevelpos) & !is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & !is.na(meanlevelUB1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind13) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	mean0 	<- meanlevelUB1[i]					# untransformed mean value
	logmean0 	<- log(mean0)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4)	# starting vector for optimization
	maxlevel0 	<- log(maxlevel[i])					# log-transformed max value
		
	hres	<- optim(par0, f.fun13, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov		<- solve(hres$hessian)				# co-variance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2]	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		FLOD		<- pnorm(LOD0, hpar[1], hpar[2])
		hres1 	<- c(13, hpar, hse, 1 - FLOD, hmeanall)		
		print(hres1)
		hest[i, 1:7] <- hres1 						# estimation results stored in output array
	}
}}



#-------------------------------------------------------------------------
# case 14: known fraction positives, LOQ, and mean level
#-------------------------------------------------------------------------

hind14 <- which(is.na(medlevel) & is.na(meanlevel) & is.na(meanlevelpos) & is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & !is.na(meanlevelUB1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind14) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	mean0all 	<- meanlevelUB1[i]					# untransformed mean value
	logmean0 	<- log(mean0all)
	par0 		<- c(logmean0 + eps, abs(logmean0 + eps) / 4)	# starting vector for optimization	
		
	hres	<- optim(par0, f.fun14, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov		<- solve(hres$hessian)				# co-variance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2]	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		FLOD		<- pnorm(LOD0, hpar[1], hpar[2])
		hres1 	<- c(14, hpar, hse, 1 - FLOD, hmeanall)		
		print(hres1)
		hest[i, 1:7] <- hres1 						# estimation results stored in output array
	}
}}

#-------------------------------------------------------------------------
# case 15: known fraction positive, mean value over positives, and LOQ
#-------------------------------------------------------------------------

hind15 <- which(is.na(medlevel) & is.na(meanlevel) & !is.na(meanlevelpos) & is.na(maxlevel) & !is.na(poscase) & !is.na(LOQ) & !is.na(nsample1) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind15) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	mean0pos 	<- meanlevelpos[i]					# untransformed mean value
	logposmean0 <- log(mean0pos)
	par0 		<- c(logposmean0 + eps, abs(logposmean0 + eps) / 4)	# starting vector for optimization	
		
	hres	<- optim(par0, f.fun15, hessian = TRUE) 			# optimization
	print(hres)
	hpar	<- hres$par
	if (hpar[1] > minlogmu) {						# acceptable model fit
		hpar[2] 	<- exp(hpar[2])					# back-transformation of sigma-parameter
		hcov		<- solve(hres$hessian)				# co-variance matrix
		hse		<- sqrt(diag(hcov))				# calculated standard errors of parameter estimates
		hcorrmusig[i] <- hcov[1, 2] / prod(hse)			# correlation between estimated mu and sig
		hse[2]	<- hpar[2] * hse[2]				# adjustment for sigma-parameter
		FLOD		<- pnorm(LOD0, hpar[1], hpar[2])
		hres1 	<- c(15, hpar, hse, 1 - FLOD, hmeanpos)		
		print(hres1)
		hest[i, c(1:6, 8)] <- hres1 						# estimation results stored in output array
	}
}}

#-------------------------------------------------------------------------
# case 16: known LOQ, and no positive cases
#-------------------------------------------------------------------------

hind16 <- which((hest[ , 1] == -1) & !is.na(LOQ) & !is.na(nsample1) & (poscase == 0) & (hest[ , 1] == -1))

if (TRUE) { for (i in hind16) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	pest0		<- .95^(1 / nsample)
	mu3		<- log(.5 * LOD[i])
	hsig1		<- c(.01, 1)
	hres		<- optimize(f.fun16, hsig1)
	sigopt	<- hres$minimum
	sigse		<- sqrt(f.calcvar(f.fun16, sigopt, .1))
	hest[i, c(1:3, 5)] <- c(16, mu3, sigopt, sigse)
}}	

#-------------------------------------------------------------------------
# case 17: known LOQ, and some positive cases
#-------------------------------------------------------------------------

hind17 <- which((hest[ , 1] == -1) & !is.na(LOQ) & !is.na(nsample1) & (poscase > 0) & (hest[ , 1] == -1))
sigopt <- mean(hest[hest[,1] == 1, 3])
sigoptse <- sqrt(var(hest[hest[,1] == 1, 3]))

if (TRUE) { for (i in hind17) {
	print(paste("row", i))
	poscase0 	<- poscase[i]
	nsample 	<- nsample1[i]
	
	LOD0 		<- log(LOD[i])						# log-transformed LOD-value
	LOQ0 		<- log(LOQ[i])						# log-transformed LOQ-value
	mu3		<- log(.5 * LOD[i])
	hres		<- optimize(f.fun17, c(.5, 2) * mu3)
	muopt		<- hres$minimum
	muse		<- sqrt(f.calcvar(f.fun17, muopt, .1))
	hest[i, c(1:5)] <- c(17, muopt, sigopt, muse, sigse)
}}	

#-------------------------------------------------------------------------
# check on all subsets being disjunct
# the matrix hintersect contains the number of rows of the excel-file that are included
# in the subsets i (row-number) and j(column-number)
# if the result is a diagonal matrix, all subsets are indeed disjunct
#-------------------------------------------------------------------------

print("rows selected")
print(sort(c(hind0, hind1, hind2, hind3, hind4, hind5, hind6, hind7, hind8, hind9, hind10, 
			hind11, hind12, hind13, hind14, hind15, hind16, hind17)))
hindsel <- list(hind0, hind1, hind2, hind3, hind4, hind5, hind6, hind7, hind8, hind9, hind10,
			hind11, hind12, hind13, hind14, hind15, hind16, hind17)
hname <- 0:(length(hindsel) - 1)
hintersect <- array(0, dim = c(length(hindsel), length(hindsel)))
dimnames(hintersect) <- list(row = hname, col = hname)
for (i in 1:length(hindsel)) { for (j in 1:length(hindsel)) { 
	hintersect[i, j] <- length(intersect(hindsel[[i]], hindsel[[j]])) }}
print(hintersect)

#-------------------------------------------------------------------------
# more on non-selected rows
#-------------------------------------------------------------------------

hindnonsel	<- which((hest[ , 1] == -1) & !(is.na(nsample1)))
print(cbind(hindnonsel,hdat[hindnonsel, c(2,3,13:14)]))

if (FALSE) {

	indmat <- cbind(!is.na(medlevel), !is.na(meanlevel), !is.na(meanlevelpos), !is.na(maxlevel), !is.na(poscase), !is.na(LOQ),
		!is.na(nsample1), !is.na(meanlevelUB1))
	indmat[91,] }

#-------------------------------------------------------------------------
# presentation of all results in one excel file, with as rows the rows of the data-excel file
# for which I could fit a log-normal distribution, and as columns the data and calculated 
# distributional characteristics
#-------------------------------------------------------------------------

if (TRUE) {

	muback	<- exp(hest[ , 2] + .5 * hest[ , 3]^2)
	sigback	<- muback * sqrt(exp(hest[ , 3]^2) - 1)

	hindall	<- which(hest[ , 1] > (-.5))
	hres		<- cbind(hest[, 1:5], muback, sigback, poscase, hest[, 6], meanlevel, hest[, 7], meanlevelpos, hest[, 8], medlevel, hest[ , 9])
		# one array with data and estimated values combined
	hres		<- hres[hindall, ]

	hres1		<- array(round(hres, 2), dim = dim(hres))
	hres1 	<- cbind(hdat1[hindall, 1:4], hres1)

#-------------------------------------------------------------------------
# estimation of missing (non-calculated) standard errors of mu and/or sigma
# 2-step procedure: first: based on same products, second: based on all products
#-------------------------------------------------------------------------

	nhres1 <- dim(hres1)[1]

	for (k in 4:5) {

		hind2 <- which((hres[ , k] < 1E-8) | is.na(hres[ , k]))

		for (i1 in 1:length(hind2)) {
			i	<- hind2[i1]
			if (!is.na(hres1[i,2])) {
				j <- which((hres1[ , 2] == hres1[i, 2]) & (hres[ , k] > 1E-8))
				hres[i , k] <- mean(hres[j, k])
			} else { 
				if (!is.na(hres1[i,3])) {
					j <- which((hres1[ , 3] == hres1[i, 3]) & (hres[ , k] > 1E-8))
					hres[i, k] <- mean(hres[j, k])
				}
			}
		}
	}

	for (k in 4:5) {
		hind2 <- which((hres[ , k] < 1E-8) | is.na(hres[ , k]))
		hind2rest <- setdiff(1:dim(hres)[1], hind2)
		hval <- log(hres[hind2rest, k] / abs(hres[hind2rest, k - 2]))
		hres[hind2, k] <- exp(mean(hval)) * abs(hres[hind2, k - 2])
	}
	hres1[ , 4 + (4:5)] <- hres[ , 4:5]

#-------------------------------------------------------------------------
# write to file
#-------------------------------------------------------------------------

	hres1		<- data.frame(hres1) # array becomes data frame
	names(hres1) <- c("type", "raw", "processed", "sample size", "method", "mu", "se", "mu_se", "se_se", "mu_orig", "sig_orig",  "poscase_dat", "poscase_est", "meanlevel_dat", "meanlevel_est", 
				"meanlevelpos_dat", "meanlevelpos_est", "medlevel_dat", "medlevel_est")
	hres1$ref 	<- refer[hindall] # add reference to data frame
	for (i in 4:17) { hres1[[i]] <- as.numeric(hres1[[i]]) }
	print(hres1[ , -dim(hres1)[2]])
	hfilename <- paste("resultstudies1", currdate, ".csv", sep = "")
	write(names(hres1), hfilename, ncol = length(hres1), sep = ";", append = FALSE)
	write(t(as.matrix(hres1)), hfilename, ncol = length(hres1), sep = ";", append = TRUE) # write results to csv-file

}


