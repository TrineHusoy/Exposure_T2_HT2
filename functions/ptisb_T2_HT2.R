# Calculation of air:blood partition coefficient
get_air_blood_partition <- function(VP, logKow, molW, waterSol) {
  #' @param VP vapour pressure [kg*m^-1*s^-2]
  #' @param lowKow Log octanol:water partition [-]
  #' @param molW molecular weight [g/mol]
  #' @param waterSol water solubility [g/L]
  gasConst <- 8.314 # kg*m^2*s^-2
  temp <- 293 # K
  henryDL <- VP * molW / (waterSol*1000 * gasConst * temp ) # unitless # divided by 1000 to convert L to m3
  logKoa <- logKow - log10(henryDL)
  
  if (VP > 4000 & henryDL > 0.1) pba <- 0.8417 / henryDL + 0.006232 * 10^logKoa
  else pba <- 0.4445 / henryDL + 0.005189 * 10^logKoa
  return(1/pba)
}

get_tissue_blood_partition <- function (logKolw=NA, fup=NA, fut, logKow=NA, ion=NA, VP=NA, molW=NA, waterSol=NA) 
{
  #' @title Helper Functions for PTISB Modelling
  #' @description This file contains helper functions used in the PTISB (Physiologically Based Kinetic) modelling workflow.
  #' @details These functions support various tasks required for PBK model development and analysis.
  #' @param logKolw Numeric value of log octanol:water partition coefficients (if NA, logKow is used as surrogate)
  #' @param fup Numeric value of unbound fractions in plasma (if NA, estimated from logKow and ionization)
  #' @param fut Numeric value of unbound fractions in tissues (if NA, estimated from fup)
  #' @param logKow Numeric value of log octanol:water partition coefficients (used if logKolw is NA)
  #' @param ion Character value indicating predominant state of ionization at pH 7.4 ("neg", "neu", "pos", "zwi", "pps", or "pos2")  
  #' @param VP Numeric value of vapour pressures [kg*m^-1*s^-2] (needed for air:blood partitioning)
  #' @param molW Numeric value of molecular weights [g/mol] (needed for air:blood partitioning)
  #' @param waterSol Numeric value of water solubilities [g/L] (needed for air:blood partitioning)
  #' @return Data frame of tissue:blood partitioning coefficients and input properties
  #' 
  #' @export
  
  if (is.na(logKolw) & is.na(logKow)) {
    stop("logKolw or logKow must be provided.")
  }
  if (is.na(fut) & is.na(ion)) {
    stop("fut (unbound fraction in tissues) must be provided if ionization is unknown.")
  }
  
  #   logKow and ion may be NA unless logKolw, fup and fut contain unknowns:
  #     (1) If logKolw is NA, it is replaced by logKow as a surrogate; 
  #     (2) NAs in fup are replaced by estimates from logKow and ionization 
  #         according to Lobell and Sivarajah (2003) or Obach (2008);
  #     (3) NAs in fut are replaced by estimates from fup according to Poulin and
  #         Theil (2000).

  Ptisb <- data.frame(Plunb = NA,
                       Phrtb = NA,
                       Psknb = NA,
                       Padpb = NA,
                       Pmusb = NA,
                       Pbonb = NA,
                       Pbrnb = NA,
                       Pthyb = NA,
                       Pgonb = NA,
                       Pkidb = NA,
                       Pstob = NA,
                       Pintb = NA,
                       Psplb = NA,
                       Ppanb = NA,
                       Plivb = NA,
                       Premb = NA,
                       Pab = NA,
                       logKolw = logKolw,
                       fup = fup,
                       fut = fut,
                       logKow = logKow,
                       ion = ion)

  
  ## -----------------------------------------------------------------------------  
  #  Matrix of tissue compositions (Poulin, from ICRP?)
  ## -----------------------------------------------------------------------------
  # cols: "Vw"  = volume fraction water 
  #       "Vn"  = volume fraction neutral lipids
  #       "Vp"  = volume fraction phospholipids            
  # rows: lun, hrt, skn, adp, mus, bon, brn, thy, gon, kid, sto, int, spl, 
  #            pan, liv, rem, pl
  
  Vw  <- c(0.78, 0.77, 0.7, 0.15, 0.74, 0.35, 0.75, 0.75, 0.77, 0.75, 0.7, 0.7, 
           0.77, 0.77, 0.75, 0.75, 0.94)
  Vn  <- c(0.0199, 0.0117, 0.0205, 0.798, 0.009, 0.0222, 0.0393, 0.01, 0.0077, 
           0.0393, 0.032, 0.032, 0.0077, 0.0077, 0.1, 0.01, 0.0032)
  Vp  <- c(0.017, 0.0141, 0.0155, 0.002, 0.01, 0.0005, 0.0532, 0.01, 0.0136, 
           0.0532, 0.015, 0.015, 0.0136, 0.0136, 0.1, 0.01, 0.002)
  
  tiscomp <- cbind(Vw, Vn, Vp)
  
  ## -----------------------------------------------------------------------------
  #  Missing information
  ## -----------------------------------------------------------------------------
  
  if (length(which(is.na(Ptisb$logKolw))) > 0)
  { 
    # Replace unknowns in logKolw by logKow (QPPR unknown!)
    logKolw.isna <- which(is.na(Ptisb$logKolw))
    Ptisb$logKolw[logKolw.isna] <- Ptisb$logKow[logKolw.isna]
  }  
  
  if (length(which(is.na(Ptisb$fup))) > 0) 
  {                       
    ## Find unknowns in fup                 
    fup.isna <- which(is.na(Ptisb$fup))
    
    if (length(which(Ptisb$ion == "pos" & Ptisb$logKow > 0.2)) > 0)
    {
      ## Split positively charged compounds into logKow > 0.2 and <= 0.2  
      Ptisb$ion[which(Ptisb$ion == "pos" & Ptisb$logKow > 0.2)] <- "pos2" 
    }
    ## Replace unknown fup by QPPR estimate from logKow and ionization
    
    f.qppr.fup <- function (x, ion)
    {
      # Lobell and Sivarajah (2003)
      # Takes as argument LogKow and the predominant state of ionization at pH 7.4
      # Returns the unbound fraction in plasma
      #TODO: vectorize for loop
      
      switch(ion,
             neg  = { 
               fup <- 1/(10^(0.3649*x + 0.4162) + 1) 
             },
             zwi  = { 
               fup <- 1/(10^(0.27*x - 0.1) + 1)        # from Obach et al. (2008)    
             }, 
             neu  = { 
               fup <- 1/(10^(0.4458*x - 0.4782) + 1) 
             },
             pos  = { 
               fup <- 1/(10^(-0.945) + 1)            
             },
             pos2 = { 
               fup <- 1/(10^(0.4628*x - 1.0971) + 1) 
             },
             pps  = { 
               fup <- 1/(10^(0.3978*x - 2.0965) + 1) 
             })
    }  
    
    for (ii in fup.isna) 
    {
      Ptisb$fup[ii] <- f.qppr.fup(Ptisb$logKow[ii], Ptisb$ion[ii])
    }
  }
  
  if (length(which(is.na(Ptisb$fut))) > 0) 
  {
    ## Find unknowns in fut 
    fut.isna <- which(is.na(Ptisb$fut))
    
    ## Replace unknown fut by estimate from fup (Poulin and Theil 2000)    
    Ptisb$fut[fut.isna] <- 1/(0.5 + 0.5/Ptisb$fup[fut.isna])
  }
  ## -----------------------------------------------------------------------------
  #  Calculation of tissue:blood partitioning coefficients
  ## -----------------------------------------------------------------------------

  f.qppr.Ptb <- function (x)
  {
    # Poulin et al. (2000, 2001, 2002)
    # Argument: index of tissue col x in Ptisb
    # Returns: vector of partitioning coefficients for tissue x (length of col x)
    
    Kolw    <- 10^(Ptisb$logKolw)	
    Vwt     <- tiscomp[x,1]   # Composition tissue
    Vnt     <- tiscomp[x,2]
    Vpt     <- tiscomp[x,3]
    Vwp     <- tiscomp[17,1]  # Composition plasma
    Vnp     <- tiscomp[17,2]
    Vpp     <- tiscomp[17,3]
    fup     <- Ptisb$fup
    fut     <- Ptisb$fut
    
    Ptb <- ((Kolw*(Vnt + 0.3*Vpt) + (Vwt + 0.7*Vpt))*fup)/
      ((Kolw*(Vnp + 0.3*Vpp) + (Vwp + 0.7*Vpp))*fut)
  }
  
  tis <- 1:16
  Ptisb[,tis] <- sapply(tis, f.qppr.Ptb)  


  if (!(is.na(VP) | is.na(molW) | is.na(waterSol))) Ptisb$Pab <- get_air_blood_partition(VP, Ptisb$logKow, molW, waterSol)


  return(Ptisb)    
}