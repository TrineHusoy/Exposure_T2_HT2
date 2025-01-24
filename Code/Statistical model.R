#--------------------------------------------------------------------------------------
# program that fits mathematical-statistical transfer model to data on 24-hour intake and
# 30-hour excretion of T2/HT2 in humans
# program consists of several parts:
# - setting constants, making selections, and defining functions
# - reading data on urinary concentration, urinary volumes, and intake data
#   experiment was done in two groups (day 1, and day 2); the day kept being incuded in the data,
#   although it was not included as a confounder in the regression model 
# - calculation of the evening intake amounts before starting the experiment (optional selection)
# - defining the fit-functions; finally a mixed effects model was fitted, to be complete also the try-out fixed
#   effects model is still included in the code; the mixed effects model contains as fixed effect parameters the 
#   population mean fabs value, and the mu and sigma of the residence time; the random effect parameter
#   is the individual heterogeneity of the fabs
# - fitting the model following several steps: (1) restriction to individuals with at least 1 positive intake and
#   excretion amount, (2) creation of dataframe, (3) defining regression formula, (4) fitting model,
#   (5) plots to validate the model fit, (6) presentation of results
#   to force realistic parameter values, the logit of the fabs-parameter, and the logarithm of the sigma-parameter of the log-normal
#   distribution of the residence time were the arguments of the fit-function
# - several sensitivity analyses are possible by changing some settings or selections; two more complex
#   sensitivity analyses are re-sampling individuals, and random instead of fixed intake amounts;
#   re-sampling individuals is not reported in the paper
# - in case of random intake amounts a Monte Carlo analysis is applied, with each run using the amounts being
#   calculated from fixed intake food amounts and randomly drawn concentration values (from log-normal distribution, Trine)
#--------------------------------------------------------------------------------------

rm(list = ls())

#--------------------------------------------------------------------------------------
# packages used
#--------------------------------------------------------------------------------------

library(openxlsx)
library(nlme)

#--------------------------------------------------------------------------------------
# selections and constants
#--------------------------------------------------------------------------------------

intervname 		<- c("AM", "PM", "even")	# names of day-intervals
logind 		<- 2					# if 1, fit is on log-transformed instead of untransformed data excretion numbers, if 2: on squareroot-values
leftcensind 	<- TRUE				# if TRUE, intake and excretion in day-interval before starting the experiment equals the one during last interval; FALSE, = 0
leftcensestind 	<- FALSE				# if TRUE, intake previous evening is estimated; if FALSE, data
distribsel		<- "lognormal"			# default distribution of residence time, alternative is "gamma"
eps 			<- 1E-6				# small number
meanintakeind	<- FALSE				# using the mean (TRUE) or median (FALSE) of the random intake values for the analysis

#-------------------------------------------------------------------
# cut-off points of the day-intervals including new ones before and after data-intervals
#-------------------------------------------------------------------

bordextra 	<- c(54, 78) # 54 means 24 hours after last measurement time-point 30 o'clock; used for sensitivity-analysis: bordextra <- altern: c(50, 70) main: c(54, 78)
bord 		<- c(-6, 6, 12, 18, 32.3, bordextra) # the value 32.3 is used to include the person that had its morning excretion time point rather late

#--------------------------------------------------------------------------------------
# functions to be used for logit (back)transformation
#--------------------------------------------------------------------------------------

f.logit <- function(x) { log(x / (1 - x)) }
f.logitinv <- function(x) { exp(x) / (1 + exp(x)) }

#-------------------------------------------------------------------
# cumulative distribution depending on the selected distribution type (distribsel)
#-------------------------------------------------------------------

if (distribsel == "lognormal") { 
	f.cumprob <- function(ht, mu1, sig1) { plnorm(ht, mu1, sig1) }}
if (distribsel == "gamma") { 
	f.cumprob <- function(ht, mu1, sig1) { pgamma(ht, mu1, rate = sig1) }}

#-------------------------------------------------------------------
# functions to extract substrings and transform times to numericals
#-------------------------------------------------------------------

f.fun1 <- function(hc) { substr(hc, 1, nchar(hc) - 2) }
f.fun2 <- function(hc) { substr(hc, 1, nchar(hc) - 8) }
f.recodetime1 <- function(t1) {
	t2 <- floor(t1 / 100)
	t3 <- (t1 / 100 - t2) * 5/3
	return(t2 + t3) }
f.recodetime2 <- function(t1) {
	t2 <- as.numeric(strsplit(t1, ":")[[1]])
	t3 <- t2[1] + t2[2] / 60
	return(t3) }

#-------------------------------------------------------------------
# urinary concentration data; unit: ng/ml
#-------------------------------------------------------------------

hfile 	<- "Urinaryconcentrations.xlsx"
hdat 		<- read.xlsx(hfile, 1)
ind1 		<- hdat[[2]][1:75]
ind1A 	<- unlist(lapply(as.list(ind1), f.fun1))
ind1conc	<- as.numeric(unlist(lapply(as.list(ind1), f.fun2)))
ind1conc	<- matrix(ind1conc, byrow = TRUE, ncol = 3)
ind2 		<- hdat[[7]][1:45]
ind2A 	<- unlist(lapply(as.list(ind2), f.fun1))
ind2conc	<- as.numeric(unlist(lapply(as.list(ind2), f.fun2)))
ind2conc	<- matrix(ind2conc, byrow = TRUE, ncol = 3)

#-------------------------------------------------------------------
# HT-2 data for days 1 and 2
#-------------------------------------------------------------------

HT21 		<- matrix(as.numeric(hdat[[4]][1:75]), byrow = TRUE, ncol = 3)
HT22 		<- matrix(as.numeric(hdat[[9]][1:45]), byrow = TRUE, ncol = 3)

#-------------------------------------------------------------------
# urinary volume data; unit: ml
#-------------------------------------------------------------------

hfile		<- "Urinaryvolumes.csv"
hdat 		<- read.csv(hfile, header = FALSE, sep = ";")
ind1 		<- hdat[[1]][1:75]
ind1A 	<- unlist(lapply(as.list(ind1), f.fun1))
ind1vol 	<- as.numeric(unlist(lapply(as.list(ind1), f.fun2)))
ind1vol 	<- matrix(ind1vol, byrow = TRUE, ncol = 3)
vol1 		<- hdat[[2]][1:75]
vol1 		<- matrix(vol1, byrow = TRUE, ncol = 3)
timeexcr1 	<- unlist(lapply(as.list(hdat[[3]]), f.recodetime2))
timeexcr1 	<- matrix(timeexcr1, byrow = TRUE, ncol = 3)
timeexcr1[ , 3] <- timeexcr1[ , 3] + 24 # 24-hour correction

ind2 		<- hdat[[5]][1:45]
ind2A 	<- unlist(lapply(as.list(ind2), f.fun1))
ind2vol 	<- as.numeric(unlist(lapply(as.list(ind2), f.fun2)))
ind2vol 	<- matrix(ind2vol, byrow = TRUE, ncol = 3)
vol2 		<- hdat[[6]][1:45]
vol2 		<- matrix(vol2, byrow = TRUE, ncol = 3)
timeexcr2 	<- unlist(lapply(as.list(hdat[[7]]), f.recodetime2))
timeexcr2 	<- matrix(timeexcr2, byrow = TRUE, ncol = 3)
timeexcr2 	<- timeexcr2[1:dim(ind2vol)[1], ]
timeexcr2[ , 3] <- timeexcr2[ , 3] + 24 # 24-hour correction

timeexcr1 <- cbind(timeexcr1[ , 3] - 24, timeexcr1)
timeexcr2 <- cbind(timeexcr2[ , 3] - 24, timeexcr2)

#-------------------------------------------------------------------
# urinary amounts data; note: array is ordered by day; division by mol-weight
#-------------------------------------------------------------------

amount1 	<- vol1 * HT21	# day 1
amount2 	<- vol2 * HT22	# day 2
amount12 	<- rbind(amount1, amount2) / 424.5

if (leftcensind) { amount12 <- cbind(amount12[ , 3], amount12)
} else { amount12 <- cbind(0, amount12) }

#-------------------------------------------------------------------
# intake data for T2 and HT2; unit: ng
# mean or median (depending on selection made above) intake values are calculated to be used in the analyses to follow
# the Monte Carlo intake values are used in the sensitivity analysis
# division by mol-weight
#-------------------------------------------------------------------

hfile		<- "IntakeT2firstday.xlsx"				# HT2
hdat		<- read.xlsx(hfile, 1)
n1 		<- dim(hdat)[1]
ID1		<- hdat$Idkode
day1		<- hdat$Day
htime1	<- hdat$Time
htime1	<- unlist(lapply(as.list(htime1), f.recodetime1))
inta1 	<- hdat[-(1:4)] / 424.5
hmean 	<- apply(log(inta1), 1, mean)
hvar 		<- apply(log(inta1), 1, var)
if (meanintakeind) { meaninta1 <- exp(hmean + .5 * hvar)
} else { meaninta1 <- exp(hmean) }
intake10	<- intake1 <- cbind(ID1, day1, htime1, meaninta1)

hfile 	<- "IntakeHT2firstday.xlsx"				# T2
hdat 		<- read.xlsx(hfile, 1)
n2 		<- dim(hdat)[1]
ID2		<- hdat$Idkode
day2		<- hdat$Day
htime2	<- hdat$Time
htime2 	<- unlist(lapply(as.list(htime2), f.recodetime1))
inta2 	<- hdat[-(1:4)] / 466.5
hmean 	<- apply(log(inta2), 1, mean)
hvar 		<- apply(log(inta2), 1, var)
if (meanintakeind) { meaninta2 <- exp(hmean + .5 * hvar)
} else { meaninta2 <- exp(hmean) }
intake20	<- intake2 <- cbind(ID2, day2, htime2, meaninta2)

#-------------------------------------------------------------------
# intake data for T2 and HT2 the other day; unit: ng
# mean or median (depending on selection made above) intake values are calculated to be used in the analyses to follow
# the Monte Carlo intake values are used in the sensitivity analysis
#-------------------------------------------------------------------

hfile		<- "IntakeT2otherday.xlsx"				# HT2
hdat		<- read.xlsx(hfile, 1)
n1 		<- dim(hdat)[1]
ID1A		<- hdat$Idkode
day1A		<- hdat$Day
htime1A	<- hdat$Time
htime1A	<- unlist(lapply(as.list(htime1A), f.recodetime1))
inta1A 	<- hdat[-(1:4)] / 424.5
hmean 	<- apply(log(inta1A), 1, mean)
hvar 		<- apply(log(inta1A), 1, var)
if (meanintakeind) { meaninta1A <- exp(hmean + .5 * hvar)
} else { meaninta1A <- exp(hmean) }
intake10A	<- intake1A <- cbind(ID1A, day1A, htime1A, meaninta1A)

hfile 	<- "IntakeHT2otherday.xlsx"				# T2
hdat 		<- read.xlsx(hfile, 1)
n2 		<- dim(hdat)[1]
ID2A		<- hdat$Idkode
day2A		<- hdat$Day
htime2A	<- hdat$Time
htime2A 	<- unlist(lapply(as.list(htime2A), f.recodetime1))
inta2A 	<- hdat[-(1:4)] / 466.5
hmean 	<- apply(log(inta2A), 1, mean)
hvar 		<- apply(log(inta2A), 1, var)
if (meanintakeind) { meaninta2A <- exp(hmean + .5 * hvar)
} else { meaninta2A <- exp(hmean) }
intake20A	<- intake2A <- cbind(ID2A, day2A, htime2A, meaninta2A)

#-------------------------------------------------------------------
# prediction of evening intake during previous day from data day-interval intake values
# whether these predicted values depends on the selections made above
#-------------------------------------------------------------------

#source("correlationsintake.R")

#--------------------------------------------------------------------------------
# unique ID values
#--------------------------------------------------------------------------------

hID 		<- unique(intake1[ , 1])
hIDA 		<- unique(intake1A[ , 1])

#-------------------------------------------------------------------
# adding evening-intake before starting experiment from the data evening-intake values
# amounts stay the same, time points minus 24 hrs
# is done for both groups of persons separately
# value may be 0, depending on selection made above
#-------------------------------------------------------------------

intake10extra <- inta1extra <- c()
for (i in 1:dim(intake10)[1]) {
	if ((intake10[i, 3] > 17) & (intake10[i, 3] < 30)) {
		intake10extra	<- rbind(intake10extra, intake10[i, ] - c(0, 0, 24, 0))
		inta1i 		<- inta1[i, ]
		hi 			<- match(round(intake10[i, 1]), hID)
		if (leftcensestind) { inta1i <- eveningratio[hi] * (inta1i + eps) - eps }
		inta1extra 		<- rbind(inta1extra, inta1i) }}
if (!leftcensind) {
	intake10extra[ , 4]	<- 0 * intake10extra[ , 4]
	inta1extra 			<- 0 * inta1extra }

inta1 	<- rbind(inta1, inta1extra)
intake10 	<- rbind(intake10, intake10extra)
hord 		<- order(intake10[ , 1], intake10[ , 3])
inta1 	<- inta1[hord, ]
intake10 	<- intake1 <- intake10[hord, ]

intake20extra <- inta2extra <- c()
for (i in 1:dim(intake20)[1]) {
	if ((intake20[i, 3] > 17) & (intake20[i, 3] < 30)) {
		intake20extra 	<- rbind(intake20extra, intake20[i, ] - c(0, 0, 24, 0))
		inta2i 		<- inta2[i, ]
		hi 			<- match(round(intake20[1]), hID)
		if (leftcensestind) { inta2i <- eveningratio[hi] * (inta2i + eps) - eps }
		inta2extra 		<- rbind(inta2extra, inta2i) }}
if (!leftcensind) {
	intake20extra[ , 4] 	<- 0 * intake20extra[ , 4] 
	inta2extra 			<- 0 * inta2extra }

inta2 	<- rbind(inta2, inta2extra)
intake20 	<- rbind(intake20, intake20extra)
hord 		<- order(intake20[ , 1], intake20[ , 3])
inta2 	<- inta2[hord, ]
intake20 	<- intake2 <- intake20[hord, ]

#-------------------------------------------------------------------
# calculation of mean or median (depending on selection made above) 
# value of intake via the mean and variance of the log-transformed values
# (since the distributiion was assumed log-normal
#-------------------------------------------------------------------

hmean 	<- apply(log(inta1 + eps), 1, mean)
hvar		<- apply(log(inta1 + eps), 1, var)
if (meanintakeind) { meaninta1 <- exp(hmean + .5 * hvar)
} else { meaninta1 <- exp(hmean) }
intake1[ , 4] <- intake10[ , 4] <- meaninta1

hmean 	<- apply(log(inta2 + eps), 1, mean)
hvar		<- apply(log(inta2 + eps), 1, var)
if (meanintakeind) { meaninta2 <- exp(hmean + .5 * hvar)
} else { meaninta2 <- exp(hmean) }
intake2[ , 4] <- intake20[ , 4] <- meaninta2

#-------------------------------------------------------------------
# unique individual IDs
#-------------------------------------------------------------------

IDuniq	<- c(ind1vol[ , 3], ind2vol[ , 3])
nIDuniq	<- length(IDuniq)
hmatch	<- match(sort(IDuniq), IDuniq)
amount12 	<- amount12[hmatch, ]
dayID		<- c(rep(1, 25), rep(2, 15))[hmatch]
IDuniq	<- sort(IDuniq)

#-------------------------------------------------------------------
# fit-function used in case of fixed effects model; NB: not used in the final model fit
# argument:
# hpar: parameter values to be estimated
# 1: fabs-value
# 2,3: mu and sigma of log-normal distribution of time-delay
# 4 (when included): other fabs-value
# result: sum-of-squares = model fit + penalty on excretions inlast interval
#-------------------------------------------------------------------

f.fitfun <- function(hpar, IDsel1, hintake1, hintake2) {
	nIDsel 	<- length(IDsel1)
	fabs1		<- fabs2 <- exp(hpar[1]) / (1 + exp(hpar[1]))			# 1st fabs-value, back-transformed
	if (length(hpar) > 3) { fabs2 <- exp(hpar[4]) / (1 + exp(hpar[4])) } 	# 2nd fabs-value, back-transformed; may be introduced because of smaller values 2nd day
	fabs 		<- c(fabs1, fabs2)[dayID[IDsel1]]					# fabs-values for all individuals
	hmu		<- hpar[2]									# mu-parameter of time-delay distribution
	hsig		<- hpar[3]									# sig-parameter of time-delay distribution

	hexcr <<- array(0, dim = c(nIDsel , 5))						# calculated excretion amounts for all individuals for all intervals

	for (i0 in 1:nIDsel ) {									# for-loop over individuals
		i1	<- IDuniq[IDsel1[i0]]
		d	<- dayID[IDsel1[i0]]
		if (d == 1) {
			i2		<- match(i1, ind1vol[ , 3])				# rows with excretion data
			time_excr	<- c(timeexcr1[i2, ], bordextra) }			# excretion time points
		if (d == 2) {
			i2		<- match(i1, ind2vol[ , 3])
			time_excr	<- c(timeexcr2[i2, ], bordextra) }

		hintake	<- which(hintake1[ , 1] == i1)				# rows with intake data
		hnintake	<- length(hintake)						# # rows
		tintake 	<- hintake1[hintake, 3]						# intake time points
		mintake 	<- hintake1[hintake, 4] + hintake2[hintake, 4]		# intake amounts

		for (i in 1:4) {									# for-loop over excretion intervals
			ind1	<- which(tintake <= time_excr[i + 1])			# intakes (foods) before excretion time point
			nind1	<- length(ind1)							# # intakes (foods)
			for (j in seq(length = nind1)) {					# for-loop over intakes (foods)
				hexcr[i0, i] <<- hexcr[i0, i] + fabs[i0] * mintake[j] * # update of total excretion in interval
					(f.cumprob(time_excr[i + 1] - tintake[j], hmu, hsig) - f.cumprob(time_excr[i] - tintake[j], hmu, hsig										))
			}
		}
		i		<- 5
		for (j in seq(length = nind1)) {						# i==5, last time-interval that includes amounts not being excreted so far
			hexcr[i0, i] <<- hexcr[i0, i] + fabs[i0] * mintake[j] * # update of total excretion in interval
				(1 - f.cumprob(time_excr[i] - tintake[j], hmu, hsig))
		}
	}
	hdiff 	<<- (hexcr[ , 1:3] - amount12[IDsel1, 2:4])			# differences between calculated and data excretion amounts
	hdiff2 	<- hdiff^2									# squared differences
	if (logind == 1) { 
		hdiff2 <- (log((hexcr[ , 1:3] + eps) / (amount12[IDsel1, 2:4] + eps)))^2 } 	# squared differences of log-transformed values
	if (logind == 2) {
		hdiff2 <- (sqrt(hexcr[ , 1:3]) - sqrt(amount12[IDsel1, 2:4]))^2 }			# squared differences of lsquareroot values

	hsom 		<- sum(hdiff2[!is.na(hdiff2)])					# sum of squared differences
	hpen 		<- 1e3 * sum(hexcr[ ,5])						# penalty on excretions in last interval
	return(hsom + hpen) }									# result is sum of squared differences + penalty

#-------------------------------------------------------------------
# fit-function in case of mixed effects model; NB: used in the final analysis
# argument:
# hpar: parameter values to be estimated
# 1: fabs-value
# 2,3: mu and sigma of log-normal distribution of residence time
# 4 (when included): other fabs-value
# result: calculated excretion amounts (may be log-transformed)
#-------------------------------------------------------------------

f.fitfunmixed <- function() {
	nIDsel 	<- length(IDsel)
	hexcr 	<<- array(0, dim = c(nIDsel , 5))				# calculated excretion amounts for all individuals for all intervals

	for (i0 in 1:nIDsel ) {								# for-loop over individuals

		hmu	<- eval(parse(text = "parmu"))[4 * i0 - 3]		# mu-parameter of residence time distribution
		hsig	<- eval(parse(text = "parsig"))[4 * i0 - 3] 		# sig-parameter of residence time distribution
		fabs1	<- eval(parse(text = "fabs"))[4 * i0 - 3]			# fabs parameter

		#print(c("pars",i0,round(c(hmu, hsig,fabs1,logind),2)))

		fabs1 <- exp(fabs1 ) / (1 + exp(fabs1 ))				# back-transformed fabs-value (because it must be in (0,1) interval
		hsig	<- exp(hsig)							# back-transformed sig-value (because it must bepositive)
		
		i1	<- IDuniq[IDsel[i0]]						# ID
		d	<- dayID[IDsel[i0]]						# day
		if (d == 1) {
			i2		<- match(i1, ind1vol[ , 3])			# rows with excretion data for day 1
			time_excr	<- c(timeexcr1[i2, ], bordextra) }		# excretion time points
		if (d == 2) {
			i2		<- match(i1, ind2vol[ , 3])			# rows with excretion data for day 2
			time_excr	<- c(timeexcr2[i2, ], bordextra) }

		hintake	<- which(intake1[ , 1] == i1)				# rows with intake data
		hnintake	<- length(hintake)					# # rows
		tintake 	<- intake1[hintake, 3]					# intake time points
		mintake 	<- intake1[hintake, 4] + intake2[hintake, 4]	# intake amounts

		for (i in 1:4) {								# for-loop over excretion intervals except last one
			ind1	<- which(tintake <= time_excr[i + 1])		# intakes (foods) before excretion time point
			nind1	<- length(ind1)						# # intakes (foods)
			for (j in seq(length = nind1)) {				# for-loop over intakes (foods)
				hexcr[i0, i] <<- hexcr[i0, i] + fabs1 * mintake[j] * # update of total excretion in interval
					(f.cumprob(time_excr[i + 1] - tintake[j], hmu, hsig) - f.cumprob(time_excr[i] - tintake[j], hmu, hsig))
			}
		}
		i	<- 5									# last excretion interval
		for (j in seq(length = nind1)) {					# for-loop over intakes (foods)
			hexcr[i0, i] <<- hexcr[i0, i] + 100 * fabs1 * mintake[j] * 	# update of total excretion in interval; factor 100 works as penalty
				(1 - f.cumprob(time_excr[i] - tintake[j], hmu, hsig))
		}
	}
	
	hexcr[ , 1:5] <<- pmax(hexcr[ , 1:5], .001)				# minimum excretion value is set to .005
	hexcr13 <<- c(t(hexcr[ , c(1:3, 5)]))					# calculated excretion amounts are ordered in the same way as data excretion amounts
	if (logind == 1) { hexcr13 <<- log(hexcr13) }				# log-transformation if selected
	if (logind == 2) { hexcr13 <<- sqrt(hexcr13) }				# log-transformation if selected


	return(hexcr13) }									# result is vector of model-calculated values

#-------------------------------------------------------------------
# fitting the model by (1) creating the data-list, (2) defining the regression model,
# and (3) applying the nlme-procedure
#-------------------------------------------------------------------

intake1[ , 4] <- meaninta1
intake2[ , 4] <- meaninta2

#-------------------------------------------------------------------
# comparisons of intake values
#-------------------------------------------------------------------

print(summary(c(unlist(inta1))))
print(summary(meaninta1))
	
#-------------------------------------------------------------------
# selection of individuals with at least 1 positive intake and excretion amount
#-------------------------------------------------------------------

posvaluesind 	<- c(FALSE, TRUE)[2]					# indicator of selecting only individuals with at least 1 positive intake and excretion amount (value TRUE)				
	
posexcr 		<- apply(amount12, 1, sum) > 0
posintake		<- (apply(table(intake1[ , 1], intake1[ , 4]), 1, sum) +
				apply(table(intake2[ , 1], intake2[ , 4]), 1, sum) > 0)
posexcrintake 	<- (1:40) * !is.na(posexcr & posintake) * posexcr * posintake
posexcrintake 	<- unname(posexcrintake[!is.na(posexcrintake)])

IDsel 		<- 1:40							# all individuals selected
if (posvaluesind) { 								# only individuals with (at least 1) positieve intake and excretion amount included
	IDsel <- IDsel0 <- (1:40)[posexcrintake] }
nIDsel 		<- length(IDsel)						# # individuals included
	
#-------------------------------------------------------------------
# creation of dataframe
#-------------------------------------------------------------------

hdat 			<- cbind(rep(IDsel, 4), c(cbind(amount12[IDsel, 2:4], .001)))
hord 			<- order(hdat[ , 1])
hdat 			<- hdat[hord, ]
hdat 			<- data.frame(hdat)
names(hdat) 	<- c("ID", "y")
if (logind == 1) { hdat$y <- log(hdat$y) }				# if fit is in log-transformed excretion amounts
if (logind == 2) { hdat$y <- sqrt(hdat$y) }				# if fit is in log-transformed excretion amounts
hdat0 		<- hdat

#-------------------------------------------------------------------
# defining the regression model including the start-vector of parameter-values
#-------------------------------------------------------------------	

if (distribsel == "lognormal") { hvect <- c(fabs = -3., parmu = 2, parsig = -1) }
if (distribsel == "gamma") { hvect <- c(fabs = log(.03), parmu = 20, parsig = log(1)) }
formals(f.fitfunmixed, envir = environment(f.fitfunmixed)) <- as.list(hvect)
hform <- as.formula(paste("y", "~", "f.fitfunmixed(", paste(names(hvect), collapse = ","), ")"))

#-------------------------------------------------------------------
# applying the nlme-function
#-------------------------------------------------------------------	

hres <- nlme(hform, data = hdat, fixed = fabs + parmu + parsig ~ 1, random = pdDiag(fabs ~ ID),
		groups = ~ID,
		start = hvect, 
		control = list(tolerance = 1E-3, minScale = 1E-3, msTol = 1E-3, pnlsTol = 1E-3, maxIter = 500),
		verbose = TRUE)
print(summary(hres))
hparoptim <-  hres$coefficients$fixed

#-------------------------------------------------------------------
# some model results: population median and mean fabs-value
#-------------------------------------------------------------------	
	
hvarfixed		<- hres$varFix[1, 1]					# variance of population fabs-value on transformed scale
hvarrand		<- sum(as.numeric(VarCorr(hres)[1:2]))		# variance of individual fabs-value on transformed scale
hhvartot		<- hvarfixed + hvarrand

medlogitfabs 	<- hparoptim[1]						# median value fabs on transformed scale
medfabspop 		<- exp(medlogitfabs) / (1 + exp(medlogitfabs))	# idem on original scale
meanfabspop		<- medfabspop + .5 * hhvartot * medfabspop * (1 - medfabspop) * (1 - 2 * medfabspop)
varmedfabspop	<- medfabspop^2 * (1 - medfabspop)^2 * hvarfixed
varmeanfabspop	<- meanfabspop^2 * (1 - meanfabspop)^2 * hvarfixed

print("population median fabs, CI bounds, and SE")
print(round(c(medfabspop + c(0, -2, 2) * sqrt(varmedfabspop), sqrt(varmedfabspop)),3))
print("population mean fabs, CI bounds, and SE")
print(round(c(meanfabspop + c(0, -2, 2) * sqrt(varmeanfabspop), sqrt(varmeanfabspop)),3))

#-------------------------------------------------------------------
# some model results: individual fabs-value
#-------------------------------------------------------------------	

logitfabsindiv 	<- hparoptim[1] + apply(hres$coefficients$random$ID, 1, sum)	# predicted individual fabs-values on transformed scale
fabsindiv		<- exp(logitfabsindiv) / (1 + exp(logitfabsindiv))
meanfabsindiv	<- mean(fabsindiv)
varfabsindiv	<- var(fabsindiv)
medfabsindiv	<- median(fabsindiv)

print("individual fabs: median, mean, SE")
print(round(c(meanfabsindiv, medfabsindiv, sqrt(varfabsindiv)),3))

#-------------------------------------------------------------------
# some model results: residence time; distinguishing between log-normal or gamma-distribution
#-------------------------------------------------------------------	
	
hsig0		<- exp(hparoptim[3])
if (distribsel == "lognormal") {
	hmu		<- exp(hparoptim[2] + .5 * hsig0^2)			# mean value of time-delay on original scale
	hsig		<- hmu * sqrt(exp(hsig0^2) - 1)			# sig value of time-delay on original scale
	hbound	<- qlnorm(c(.5, .025, .975), hparoptim[2], hsig0) }	# median and 95% boundary values of time-delay on original scale
if (distribsel == "gamma") {
	hmu		<- hparoptim[2] / hsig0					# mean value of time-delay on original scale
	hsig		<- sqrt(hmu / hsig0^2)					# sig value of time-delay on original scale
	hbound	<- qgamma(c(.5, .025, .975), hparoptim[2], rate = hsig0) }	# 95% boundary values of time-delay on original scale

musigtim 	<- c("mean" = hmu, "median" = hbound[1], "2.5%" = hbound[2], "97.5%" = hbound[3], "sig" = hsig)
print("residence time characteristics")
print(round(musigtim, 2))
	
plot(hdat$y, hexcr13, type = "p", xlab = "data", ylab = "est")	# plot of predicted excretion amounts regressed on data amounts
hresid <- hres$residuals[ , 1]
plot(hdat$y, hresid, type = "p")						# plot of residuals regressed on data amounts, to check for systematic effects
qqnorm(hresid)									# qq-plot to check whether the residuals are indeed normally distributed

#-------------------------------------------------------------------
# doing uncertainty analysis by re-sampling individuals (1),
# and/or by taking random instead of mean intake values (2)
#-------------------------------------------------------------------

logind	<- (0:2)[3]								# bootstrap with dependent variable being the untransformed excretion amounts
posvaluesind <- c(FALSE, TRUE)[2]						# again, only including individuals with at least 1 positieve intake and excretion amount value
MCtype0	<- c(2)								# only 1 uncertainty analysis, i.e. random instead of mean intake values
nMC0 		<- 20									# # Monte Carlo runs

if (TRUE) { for (MCtype in MCtype0) {					# for-loop over all uncertainty analyses selected
	if (MCtype == 1) { 							# i.c.o. re-sampling individuals
		nMC <- min(100, nMC0)
		intake1[ , 4] <- meaninta1
		intake2[ , 4] <- meaninta2 }
	if (MCtype == 2) { 							# i.c.o. random intake amount values
		nMC <- min(1000, nMC0)
		IDsel <- IDsel0 }

	for (i2 in 1:1) {
		print(paste("MCtype", MCtype, "i2", i2))
		if (i2 == 1) { parstrt 	<- c(log(.05/(1-.05)), 1.5, 1) } #, log(.05/(1-.05))) } # initial value of function parameters

		#-------------------------------------------------------------------
		# initialization of arrays to store results
		#-------------------------------------------------------------------

		missingexcr	<- rep(0, nMC)						# proportion of excretion after stopping experiment (missing)
		hresultsall	<- vector("list", length = nMC)			# list including all model fits

		nonMC <- c(2,39,43)
		
		for (MC in (1:nMC)[-nonMC]) { 

			#-------------------------------------------------------------------
			# re-sampling individuals
			#-------------------------------------------------------------------
		  
			print(c("MC", MC))

			if (MCtype == 1) {
				IDsel <- sample(IDsel0, length(IDsel0), replace = TRUE)	
				hdat <- c()
				for (j in IDsel) {
					hdat <- rbind(hdat, cbind(j, hdat0[hdat0$ID == j, 2])) }
				hdat <- data.frame(hdat)
				names(hdat) <- names(hdat0)
			}
			#-------------------------------------------------------------------
			# random intake amounts
			#-------------------------------------------------------------------			

			if (MCtype == 2) { 
				hdat <- hdat0
				intake1[ , 4] <- inta1[ , MC]
				intake2[ , 4] <- inta2[ , MC]
			}

			#-------------------------------------------------------------------
			# fit the model; note: may not converge, so included in try-function
			#-------------------------------------------------------------------			

			hres <- try(nlme(hform, data = hdat, fixed = fabs + parmu + parsig ~ 1, random = pdDiag(fabs ~ ID),
					groups = ~ID,
					start = hvect, 
					control = list(tolerance = .2, minScale = .2, msTol = .2, pnlsTol = .2, maxIter = 50)))
			if (length(hres) == 1) { print(hres[1])							# model fit failure
			} else {
				hresultsall[[MC]] <- hres
				print(summary(hres))
	
				missingexcr[MC]	<- sum(hexcr[ , 4:5]) / sum(hexcr)				# proportion of excretion being missing
			}
		}

		hind	<- which(unlist(lapply(hresultsall, length) > 0))

		hfixedmat	<- array(0, dim = c(nMC, length(hvect)))						# calculated fixed-effect parameter values

		hvarfixed	<- hvarrand <- hhvartot <- medlogitfabs <- medfabspop <-
			meanfabspop <- varmedfabspop <- varmeanfabspop <- meanfabsindiv <-
			varfabsindiv <- medfabsindiv <- 
			htimemean	<- htimesig <- htimemed <- rep(0, nMC)					# calculated distributional characteristics of residence time (log-normally distributed) 

		for (MC in hind) {
			hres	<- hresultsall[[MC]]

			hfixedmat[MC, ] 	<- hres$coefficients$fixed
	
			hvarfixed[MC]	<- hres$varFix[1, 1]							# variance of population fabs-value on transformed scale
			hvarrand[MC]	<- sum(as.numeric(VarCorr(hres)[1:2]))				# variance of individual fabs-value on transformed scale
			hhvartot[MC]	<- hvarfixed[MC] + hvarrand[MC]

			medlogitfabs[MC] 	<- hfixedmat[MC, 1]							# median value fabs on transformed scale
			medfabspop[MC] 	<- exp(medlogitfabs[MC]) / (1 + exp(medlogitfabs[MC]))	# idem on original scale
			meanfabspop[MC]	<- medfabspop[MC] + .5 * hhvartot[MC] * medfabspop[MC] * (1 - medfabspop[MC]) * (1 - 2 * medfabspop[MC])
			varmedfabspop[MC]	<- medfabspop[MC]^2 * (1 - medfabspop[MC])^2 * hvarfixed[MC]
			varmeanfabspop[MC] <- meanfabspop[MC]^2 * (1 - meanfabspop[MC])^2 * hvarfixed[MC]

			#-------------------------------------------------------------------
			# some model results: individual fabs-value
			#-------------------------------------------------------------------	

			logitfabsindiv 	<- hfixedmat[MC, 1] + apply(hres$coefficients$random$ID, 1, sum)	# predicted individual fabs-values on transformed scale
			fabsindiv		<- exp(logitfabsindiv) / (1 + exp(logitfabsindiv))
			meanfabsindiv[MC]	<- mean(fabsindiv)
			varfabsindiv[MC]	<- var(fabsindiv)
			medfabsindiv[MC]	<- median(fabsindiv)
			htimemed[MC]	<- exp(hfixedmat[MC, 2]) 						# back-transformed mean residence time
			hsig			<- exp(hfixedmat[MC, 3])
			htimemean[MC]	<- exp(hfixedmat[MC, 2] + .5 * hsig^2) 				# back-transformed mean residence time
			htimesig[MC]	<- htimemean[MC] * sqrt(exp(hsig^2) - 1) 				# sig value of time-delay on original scale
		}

		print(c("medfabspop", round(c(mean(medfabspop[hind]), sqrt(var(medfabspop[hind]))), 3)))
		print(c("meanfabspop", round(c(mean(meanfabspop[hind]), sqrt(var(meanfabspop[hind]))), 3)))
		print(c("medfabsindiv", round(c(mean(medfabsindiv[hind]), sqrt(var(medfabsindiv[hind]))), 3)))
		print(c("meanfabsindiv", round(c(mean(meanfabsindiv[hind]), sqrt(var(meanfabsindiv[hind]))), 3)))
		print(c("htimemed", round(c(mean(htimemed[hind]), sqrt(var(htimemed[hind]))), 3)))
		print(c("htimemean", round(c(mean(htimemean[hind]), sqrt(var(htimemean[hind]))), 3)))	

	}		
}}	



