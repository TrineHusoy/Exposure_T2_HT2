f.physB <- function(xx) 
{
  ## -----------------------------------------------------------------------------
  # Physiologically based model parameters and their between-subject variability 
  # predicted from anthropometric values according to Bosgra et al., 2012. DOI: 10.3109/10408444.2012.709225
  # (as documented in "Methods - Population physiology model")
  #
  # Arguments: xx: anthropometry with columns "Sex", "Age", "Weight" and "Height"
  #            dref (refman.ar): array with dimensions [Parameter, Age, Sex]
  #                             Age in classes 1:6, sd in 7
  
  ## -----------------------------------------------------------------------------
  
  # Assert correct input
  if (!is.numeric(xx$Sex) || !(xx$Sex %in% c(1, 2))) stop("sex must be 1 (male) or 2 (female).")
  if (!is.numeric(xx$Age) || xx$Age < 0) stop("age must be a non-negative number.")
  if (!is.numeric(xx$Weight) || xx$Weight <= 0) stop("weight must be a positive number.")
  if (!is.numeric(xx$Height) || xx$Height <= 0) stop("height must be a positive number.")
  
  # Read reference data
  refman_male <- read.csv("Data/refman_male.csv", row.names = 1)
  refman_female <- read.csv("Data/refman_female.csv", row.names = 1)

  # Merge into a 3D array: [Parameter, Age, Sex]
  dref <- array(NA, dim = c(nrow(refman_male), ncol(refman_male), 2))
  dref[,,1] <- as.matrix(refman_male)
  dref[,,2] <- as.matrix(refman_female)
  dimnames(dref) <- list(rownames(refman_male), 
                         lifenames <- c('1','2','3','4','5', '6', '7'), #these indicate ages
                         gendernames <- c('1', '2')) #these indicate gender, c("Male", "Female"))
  
  yy <- list(Sex = xx$Sex, Age = xx$Age) 
  
  # Set sex and age in xx to classes, in order to subscribe refman array
  
  xx[which(yy$Age < 12), 1] <- 1     # no distinction m/f below 12 y
  xx[which(yy$Age < 0.5), 2] <- 1 
  xx[which(yy$Age >= 0.5 & yy$Age < 2.5), 2] <- 2
  xx[which(yy$Age >= 2.5 & yy$Age < 7), 2] <- 3
  xx[which(yy$Age >= 7 & yy$Age < 12), 2] <- 4
  xx[which(yy$Age >= 12 & yy$Age < 18), 2] <- 5
  xx[which(yy$Age >= 18), 2] <- 6
  
  Lxx <- length(xx$Age)
  
  # parnames <- c("BW", "H", "BMI", "Qc", "Vlun", "Vhrt", "Vskn", "Vadp", "Vmus",
  #           "Vbon", "Vbrn", "Vthy", "Vgon", "Vkid", "Vsto", "Vint", "Vspl", 
  #           "Vpan", "Vliv", "Vrem", "Vb", "QPhrt", "QPskn", "QPadp", "QPmus", 
  #           "QPbon", "QPbrn", "QPthy", "QPgon", "QPkid", "QPsto", "QPint",  
  #           "QPspl", "QPpan", "QPliv", "QPrem", "Vtid", "Qp", "Vfrc", "Sex", 
  #           "Age", "Qc.a")  
  
  physB <- data.frame(BW   = xx$Weight,       #  1
                      H    = xx$Height,       #  2
                      BMI  = xx$Weight/((0.01*xx$Height)^2),   #  3
                      Qc   = rep(NA, Lxx),    #  4
                      Vlun = rep(NA, Lxx),    #  5
                      Vhrt = rep(NA, Lxx),    #  6
                      Vskn = rep(NA, Lxx),    #  7
                      Vadp = rep(NA, Lxx),    #  8
                      Vmus = rep(NA, Lxx),    #  9
                      Vbon = rep(NA, Lxx),    # 10
                      Vbrn = rep(NA, Lxx),    # 11
                      Vthy = rep(NA, Lxx),    # 12
                      Vgon = rep(NA, Lxx),    # 13
                      Vkid = rep(NA, Lxx),    # 14
                      Vsto = rep(NA, Lxx),    # 15
                      Vint = rep(NA, Lxx),    # 16
                      Vspl = rep(NA, Lxx),    # 17
                      Vpan = rep(NA, Lxx),    # 18
                      Vliv = rep(NA, Lxx),    # 19
                      Vrem = rep(NA, Lxx),    # 20
                      Vb   = rep(NA, Lxx),    # 21
                      Qhrt = rep(NA, Lxx),    # 22
                      Qskn = rep(NA, Lxx),    # 23
                      Qadp = rep(NA, Lxx),    # 24
                      Qmus = rep(NA, Lxx),    # 25
                      Qbon = rep(NA, Lxx),    # 26
                      Qbrn = rep(NA, Lxx),    # 27
                      Qthy = rep(NA, Lxx),    # 28
                      Qgon = rep(NA, Lxx),    # 29
                      Qkid = rep(NA, Lxx),    # 30
                      Qsto = rep(NA, Lxx),    # 31
                      Qint = rep(NA, Lxx),    # 32
                      Qspl = rep(NA, Lxx),    # 33
                      Qpan = rep(NA, Lxx),    # 34
                      Qliv = rep(NA, Lxx),    # 35
                      Qrem = rep(NA, Lxx),    # 36
                      Vtid = rep(NA, Lxx),    # 37
                      Vfrc = rep(NA, Lxx),    # 38
                      fp   = rep(NA, Lxx),    # 39
                      Qp   = rep(NA, Lxx),    # 40
                      Valv = rep(NA, Lxx),    # 41
                      Sex  = yy$Sex,          # 42
                      Age  = yy$Age,          # 43
                      Qc.a = rep(NA, Lxx),    # 44
                      sigVorg = rep(NA, Lxx), # 45
                      BSA = rep(NA, Lxx) )    # 46
  ## -----------------------------------------------------------------------------
  # 1. Organ volumes
  ## -----------------------------------------------------------------------------
  
  ## Organ volumes allometrically scaled to H
  
  Ht <- xx[[4]]
  
  ## lun 
  physB$Vlun <- exp(2.0995*log(Ht) - 11.761)
  
  ## hrt
  physB$Vhrt <- exp(2.1291*log(Ht) - 12.307)
  
  ## bon
  physB$Vbon <- exp(2.6683*log(Ht) - 12.219)
  
  ## kid
  physB$Vkid <- exp(1.9277*log(Ht) - 11.183)
  
  ## sto
  physB$Vsto <- exp(2.4542*log(Ht) - 14.568)
  
  ## int 
  physB$Vint <- exp(2.4687*log(Ht) - 12.72)
  
  ## spl
  physB$Vspl <- exp(2.1581*log(Ht) - 13.061)
  
  ## pan
  physB$Vpan <- exp(2.4323*log(Ht) - 14.632)
  
  ## liv
  physB$Vliv <- exp(1.9822*log(Ht) - 9.8071)
  
  ## rem
  physB$Vrem <- exp(2.3905*log(Ht) - 12.044)
  
  ## Organs dependent on age, gender
  
  ## Brain: two-phase relationship with age, sex as background covariate
  
  MM <- which(yy$Sex == 1)
  FF <- which(yy$Sex == 2) 
  
  physB$Vbrn[MM] <- 0.405*(3.68 - 2.68*exp(-yy$Age[MM]/0.89))* 
    (exp(-yy$Age[MM]/629))
  physB$Vbrn[FF] <- 0.373*(3.68 - 2.68*exp(-yy$Age[FF]/0.89))* 
    (exp(-yy$Age[FF]/629))                       
  
  ## thy: two-phase relationship with age
  
  physB$Vthy <- 0.014*(7.1 - 6.1*exp(-yy$Age/11.9))*
    (0.14 + 0.86*exp(-yy$Age/10.3))
  
  ## gon: sex-dependent increase with age
  
  physB$Vgon[MM] <- 0.0033 + 0.053*(1-(exp(-(yy$Age[MM]/17.5)^5.4)))
  physB$Vgon[FF] <- 0.0033 + 0.090*(1-(exp(-(yy$Age[FF]/16.8)^6.7)))
  
  ## Organs dependent on BSA (Dubois and Dubois)
  
  BWt   <- physB$BW
  physB$BSA <- 0.007184*(BWt^0.425)*(Ht^0.725) 
  
  ## skn: exponentially related to BSA
  physB$Vskn <- exp(1.6375*physB$BSA - 1.9262)
  
  ## Vb: allometric scaling to BSA
  physB$Vb[MM] <- 3.3272*physB$BSA[MM] - 0.8105
  physB$Vb[FF] <- 2.6614*physB$BSA[FF] - 0.4636
  
  ## -----------------------------------------------------------------------------
  # Residual variation
  
  Ch <- which(yy$Age < 18)
  Ad <- which(yy$Age >= 18)
  
  # mean and sd of log of BW | H
  mu.BW  <- rep(NA, Lxx)
  mu.BW[Ch] <- -8.122 + 2.356*log(physB$H[Ch])
  mu.BW[Ad] <- -4.098 + 1.630*log(physB$H[Ad])
  s.BW <- 0.1547
  
  lnBW.ii <- log(BWt)
  
  # vector of correlations with BW | H
  
  rho  <- c(Vlun = 0.707,  # Vlun #  5
            Vhrt = 0.707,  # Vhrt #  6
            Vskn = 0,      # Vskn #  7
            Vadp = 0,      # Vadp #  8
            Vmus = 0.224,  # Vmus #  9
            Vbon = 0.224,   # Vbon # 10
            Vbrn = 0,      # Vbrn # 11
            Vthy = 0.224,   # Vthy # 12
            Vgon = 0.224,   # Vgon # 13
            Vkid = 0.224,   # Vkid # 14
            Vsto = 0.707,  # Vsto # 15
            Vint = 0.224,   # Vint # 16
            Vspl = 0.224,   # Vspl # 17
            Vpan = 0.224,   # Vpan # 18
            Vliv = 0.224,   # Vliv # 19
            Vrem = 0.224,   # Vrem # 20
            Vb = 0,      # Vb   # 21
            Qhrt = 0,      # Qhrt # 22
            Qskn = 0,      # Qskn # 23
            Qadp = 0,      # Qadp # 24
            Qmus = 0,      # Qmus # 25
            Qbon = 0,      # Qbon # 26
            Qbrn = 0,      # Qbrn # 27
            Qthy = 0,      # Qthy # 28
            Qgon = 0,      # Qgon # 29
            Qkid = 0,      # Qkid # 30
            Qsto = 0,      # Qsto # 31
            Qint = 0,      # Qint # 32
            Qspl = 0,      # Qspl # 33
            Qpan = 0,      # Qpan # 34
            Qliv = 0,      # Qliv # 35
            Qrem = 0,      # Qrem # 36
            Vtid = 0,      # Vtid # 37
            vfrc = -0.707, # Vfrc # 38
            fp = 0,      # fp   # 39
            Qp = 0,      # Qp   # 40
            Valv = 0 )     # Valv # 41
  
  f.var.physB <- function (x) 
  {
    # means and standard deviations of log of Vorg
    
    mu.org <- log(physB[[x]])
    s.org  <- dref[ cbind(rep(x, Lxx), rep(7, Lxx), xx[[1]]) ]
    
    rhox <- rho[[x]]
    
    # mean and sd [ log(Vorg) | log(BW) = log(BW.ii) ] 
    
    mu.org.BW <- mu.org + (rhox*s.org/s.BW)*(lnBW.ii - mu.BW)
    s.org.BW <- sqrt((1 - rhox^2)*s.org^2) 
    
    # random variation
    
    Vorg.ii <- exp(rnorm(Lxx, mu.org.BW, s.org.BW))
    
    return(Vorg.ii)
  }
  
  #Vsds <- c(5:7, 10:21)  # vector of Vorgs in physB to add variation
  Vsds <- c('Vlun','Vhrt','Vskn','Vbon','Vbrn','Vthy','Vgon','Vkid','Vsto','Vint','Vspl','Vpan','Vliv','Vrem','Vb')
  physB[,Vsds] <- sapply(Vsds, f.var.physB)
  #physB$Vsds <- sapply(Vsds, f.var.physB)
  
  ## -----------------------------------------------------------------------------
  
  ## Muscle: included in allometric scaling
  
  #  ## Muscle: Fraction of fat free mass (FFM) (Price et al.) 
  #  #    FFM = sigVorg + Vmus 
  #  #    Vmus = Fmus*FFM (where Fmus: muscle as fraction of FFM, age-dependent)
  #  #         = sigVorg*Fmus/(1 - Fmus)
  #  
  # Sum of all compartments minus muscle and adipose tissue. Tissues not included in the
  # model comprise around 9% of BW, notably lymph and GI lumen content.
  
  physB$sigVorg <- physB$Vlun + physB$Vhrt + physB$Vskn + physB$Vbon + 
    physB$Vbrn + physB$Vthy + physB$Vgon + physB$Vkid + 
    physB$Vsto + physB$Vint + physB$Vspl + physB$Vpan + 
    physB$Vliv + physB$Vrem + physB$Vb + 0.09*BWt
  
  #    FatMus <- BWt - physB$sigVorg     
  
  #      AM <- which(yy$Age >= 0 & yy$Age < 2.5 & xx$Sex == 1)
  #      BM <- which(yy$Age >= 2.5 & yy$Age < 7 & xx$Sex == 1)
  #      CM <- which(yy$Age >= 7 & yy$Age < 12 & xx$Sex == 1)
  #      DM <- which(yy$Age >= 12 & yy$Age < 18 & xx$Sex == 1)
  #      EM <- which(yy$Age >= 18 & yy$Sex == 1)
  
  #      AF <- which(yy$Age >= 0 & yy$Age < 2.5 & xx$Sex == 2)
  #      BF <- which(yy$Age >= 2.5 & yy$Age < 7 & xx$Sex == 2)
  #      CF <- which(yy$Age >= 7 & yy$Age < 12 & xx$Sex == 2)
  #      DF <- which(yy$Age >= 12 & yy$Age < 18 & xx$Sex == 2)
  #      EF <- which(yy$Age >= 18 & yy$Sex == 2)
  
  
  Fmus <- rep(NA, Lxx)
  Fmus[which(yy$Sex == 1 & yy$Age < 18)] <- 0.3 + 
    0.0133*yy$Age[which(yy$Sex == 1 & yy$Age < 18)]
  Fmus[which(yy$Sex == 2 & yy$Age < 18)] <- 0.3 + 
    0.01*yy$Age[which(yy$Sex == 2 & yy$Age < 18)]                          
  Fmus[which(yy$Sex == 1 & yy$Age >= 18)] <- 0.54
  Fmus[which(yy$Sex == 2 & yy$Age >= 18)] <- 0.489
  
  physB$Vmus <- physB$sigVorg*Fmus/(1 - Fmus)
  
  # random variation muscle given Age 
  
  ChM <- which(yy$Sex == 1 & yy$Age < 18)
  ChF <- which(yy$Sex == 2 & yy$Age < 18)
  AdM <- which(yy$Sex == 1 & yy$Age >= 18)
  AdF <- which(yy$Sex == 2 & yy$Age >= 18)
  
  mu.BW[ChM] <- 2.437 + 0.111*yy$Age[ChM]
  mu.BW[ChF] <- 2.429 + 0.108*yy$Age[ChF]
  mu.BW[AdM] <- 4.368 + 0.00027*yy$Age[AdM]
  mu.BW[AdF] <- 4.118 + 0.00108*yy$Age[AdF]
  s.BW <- 0.166
  
  physB$Vmus  <- f.var.physB(9)
  
  ## Adipose tissue: BW minus sum of all other compartments, density 0.9 kg/L.
  
  physB$Vadp <- (physB$BW - (physB$Vmus + physB$sigVorg))/0.9
  
  ## -----------------------------------------------------------------------------
  # 2. Circulation
  ## -----------------------------------------------------------------------------
  
  ## Blood flows: Qorg (L/h) = QPorg*Vorg (L/h/L org * L)
  
  #QPs <- 22:36 
  QPs <- c('Qhrt','Qskn', 'Qadp', 'Qmus', 'Qbon', 'Qbrn', 'Qthy', 'Qgon', 'Qkid', 'Qsto','Qint','Qspl','Qpan','Qliv','Qrem') 
  
  # index of perfusion rates in dref
  # QPs - 16: index of organ volumes in physB
  
  f.QP <- function (x) 
  {
    # Argument x: QPs
    # returns: organ blood flows
    QPorg <- dref[ cbind(rep(x, Lxx), xx[[2]], xx[[1]]) ]
    # Correction for higher Qc/BW children (0 - 17 y)
    # Function derived from curve fitting in PROAST to correction factors  
    # obtained by f.corr.CO
    fcor.M <- rep(1, Lxx)
    fcor.M[ChM] <- 1.44*(3.26 - 2.26*exp(-(yy$Age[ChM]/13.8)^1))* 
      (exp(-(yy$Age[ChM]/14.7)^1))
    fcor.F <- rep(1, Lxx)
    fcor.F[ChF] <- 1.32*(3.26 - 2.26*exp(-(yy$Age[ChF]/13.8)^1))* 
      (exp(-(yy$Age[ChF]/14.7)^1))
    
    QPorg[ChM] <- QPorg[ChM]*fcor.M[ChM]
    QPorg[ChF] <- QPorg[ChF]*fcor.F[ChF]
    #Vorg  <- physB[[x - 16]]
    Vorg <- physB[(sub('.','V',x))] #replace 'Q' with 'V'
    Qorg  <- QPorg*Vorg
    return(Qorg)
  }
  
  physB[,QPs] <- sapply(QPs, f.QP)
  
  # Perfusion per kg organ, not per L. Density adipose tissue: 0.9 kg/L
  physB$Qadp <- 0.9*physB$Qadp
  
  ## Cardiac output: defined as sum of all regional blood flows (comprising 95%
  #  of Qc), allometric scaling as check (Qc.a)
  
  physB$Qc <- ( physB$Qhrt + physB$Qskn + physB$Qadp + physB$Qmus + 
                    physB$Qbon + physB$Qbrn + physB$Qthy + physB$Qgon + 
                    physB$Qkid + physB$Qsto + physB$Qint + physB$Qspl + 
                    physB$Qpan + physB$Qliv + physB$Qrem )/0.95  
  
  physB$Qc.a <- exp(1.9894*log(Ht) - 4.22)
  
  ## -----------------------------------------------------------------------------
  # 3. Respiration
  ## -----------------------------------------------------------------------------
  
  ## -----------------------------------------------------------------------------
  # Allometric scaling to BW
  # i = individual, m = mean)
  #  
  #  BWt <- xx[[3]]
  #  BW.m <- dref[ cbind(rep(1, Lxx), xx[[2]], xx[[1]]) ]
  #  
  #  f.alloBW <- function (x)
  #    {
  #    # Argument: index of reference parameters in dref,
  #    # returns: target parameter allometrically scaled to BW
  #      BW.allo <- (BWt/BW.m)^pp
  #      par.m <- dref[ cbind(rep(x, Lxx), xx[[2]], xx[[1]]) ]
  #      par.i <- BW.allo*par.m
  #      return(par.i)
  #   }
  ## -----------------------------------------------------------------------------
  
  # Vtid # 37
  # Vfrc # 38
  # fp   # 39
  # Qp   # 40
  # Valv # 41
  
  #(1)# Tidal volume (Vtid, L): allometrically scaled to H
  
  physB$Vtid[MM] <- exp(2.27*log(Ht[MM]) - 12.2)
  physB$Vtid[FF] <- exp(1.95*log(Ht[FF]) - 10.7)
  
  
  #(2)# Ventilation rate (Qp, L/h): proportional to Qc
  
  physB$Qp[MM] <-  0.96*physB$Qc[MM] + 129
  physB$Qp[FF] <-  0.62*physB$Qc[FF] + 179 
  
  #(3)# Breathing frequency (fp, /h): given by Qp/Vtid
  
  physB$fp <- physB$Qp/physB$Vtid
  
  #(4)# Functional residual capacity (Vfrc, L)
  
  physB$Vfrc <- exp(3.1474*log(Ht) - 15.089)
  # TODO: residuals, corr. BW
  
  #(5)# Alveolar volume (Valv, L): given by Vfrc + 0.5*Vtid
  
  physB$Valv <- physB$Vfrc + 0.5*physB$Vtid
  
  ## ----------------------------------------------------------------------------- 
  # Output
  ## -----------------------------------------------------------------------------
  
  return(physB)
}

