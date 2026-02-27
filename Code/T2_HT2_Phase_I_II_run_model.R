#  Running the Physiologically Based Kinetic Model for T-2, HT-2 and their metabolites 
#  with a representative exposure for 10 adult subjects
#  
#  written by J. Westerhout, PhD (joost.westerhout@rivm.nl)
#  National Institute for Public Health and the Environment (RIVM, NL)
#  Based on the PhysB model from Bosgra et al., 2012. DOI: 10.3109/10408444.2012.709225

setwd("/home/westerj/T-2 and HT-2 supplemental material")
rm(list=ls()) # to clear out the global environment
library(deSolve)
library(ggplot2)
library(RColorBrewer)
library(paletteer)
library(tidyverse)

# Load function
source("functions/activity.R")
source("functions/activityreg.R")
source("functions/parms_pbk_T2_HT2_Phase_I_II.R")
source('functions/pbk_T2_HT2_Phase_I_II.R')
source("functions/physB.R")
source("functions/ptisb_T2_HT2.R")
source("functions/utilities.R")

# Set exposure data ----
Average_exposure <- data.frame("ID"=c(1,1,1,1,1,1),
                               "Time_PBK"=c(0,0,4,4,10,10),
                               "mean"=c(463,660,463,660,463,660),
                               "toxin"=c("T2","HT2","T2","HT2","T2","HT2"))

# Simulate the average exposure
Total_exposure <- Average_exposure

IDs <- unique(Total_exposure$ID)

# Set urine data ----
Average_urine <- data.frame("ID"=c(1),
                            "Time_PBK"=c(23.14),
                            "Cumulative_amount_HT2"=c(286.4))

Urinedata <- Average_urine

# Simulation settings ----
sim_time <- 24 # Simulation duration [h]
n_step <- 10 # Number of time steps [/h]
times <- seq(0, sim_time, by = 1/n_step)

n_samples <- 10 # length(IDs) # increase if multiple monte carlo samples need to be drawn. 
IDs <- seq(1:n_samples)
include_lung=TRUE
events=NULL

# Kinetic parameters T-2 & HT-2 and their metabolites ----
# Absorption rate constants
# We assume that the rate constants are the same for the stomach and intestine compartment
# We assume that the rate constants are the same for the parent compounds and their metabolites
ka_T2 <- 0.5 # increase ka to correct for conversion into hydroxy
ka_HT2 <- 0.5 # 0.075 # increase ka to correct for conversion into hydroxy

# Metabolic conversion rate constants
# T-2 -> metabolites
# stomach
km_sto_T2_HT2 <- 0 # h-1
km_sto_T2_T2_Phase_I <- 0 # 0 # h-1
km_sto_T2_T2_Phase_II <- 0 # h-1
# intestinal tissue
km_int_T2_HT2 <- 0.3 # h-1
km_int_T2_T2_Phase_I <- 0.1 # h-1
km_int_T2_T2_Phase_II <- 0 # h-1
# liver
CLliv_T2_HT2 <- 150 # L/h 
CLliv_T2_T2_Phase_I <- 300  #L/h 
CLliv_T2_T2_Phase_II <- 0.1 # L/h
# kidney
CLkid_T2 <- 7.5 # L/h 

# HT-2 -> metabolites
# stomach
km_sto_HT2_HT2_Phase_I <- 0.2 # h-1
km_sto_HT2_HT2_Phase_II <- 0 # h-1
# intestinal tissue
km_int_HT2_HT2_Phase_I <- 0 # h-1
km_int_HT2_HT2_Phase_II <- 0 # h-1
# liver
CLliv_HT2_HT2_Phase_I <- 100 # L/h
CLliv_HT2_HT2_Phase_II <- 18 # L/h
# kidney
CLkid_HT2 <- 7.5 # L/h

# T-2-Phase I -> excretion
# kidney
CLkid_T2_Phase_I <- 7.5 # L/h

# T-2-Phase II -> T-2-Phase I or HT-2-Phase II
# stomach
km_sto_T2_Phase_II_HT2_Phase_II <- 0 # h-1
# intestinal tissue
km_int_T2_Phase_II_HT2_Phase_II <- 0 # h-1
# liver
CLliv_T2_Phase_II_HT2_Phase_II <- 0 # L/h
# kidney
CLkid_T2_Phase_II <- 7.5 # L/h

# HT-2-Phase I -> excretion
# kidney
CLkid_HT2_Phase_I <- 7.5 # L/h

# HT-2-Phase II -> excretion
# kidney
CLkid_HT2_Phase_II <- 7.5 # L/h

# Simulating population ----
set.seed(1234)
sex <- sample(1:2,n_samples,replace=T)  # {male: 1; female: 2}
age <- sample(18:79, n_samples, replace=T) # years
weight <- sample (50:100, n_samples, replace=T) # kg
height <- sample (150:200, n_samples, replace=T) # cm

physiology_df <- data.frame(Sex=sex,
                            Age=age,
                            Weight=weight,
                            Height=height)

# Define phys-chem properties ----
physchem_T2 <- data.frame(logKolw = c(NA),
                          fup = c(NA),
                          fut = c(NA),
                          logKow = c(0.9), # https://pubchem.ncbi.nlm.nih.gov/compound/T-2-Toxin#section=Computed-Properties
                          ion = c("neu"), # charge at pH 7.4
                          VP=3.1E-11,
                          molW=466.5,
                          waterSol=NA) 

physchem_HT2 <- data.frame(logKolw = c(NA),
                           fup = c(NA),
                           fut = c(NA),
                           logKow = c(0.4), # https://pubchem.ncbi.nlm.nih.gov/compound/10093830#section=Computed-Properties
                           ion = c("neu"), # charge at pH 7.4
                           VP=NA,
                           molW=424.21,
                           waterSol=NA) 

physchem_T2_Phase_II <- data.frame(logKolw = c(NA),
                                   fup = c(NA),
                                   fut = c(NA),
                                   logKow = c(-1.1), # Unknown, but expected to be lower than parent compound. Used a delta of 2 based on https://www.sciencedirect.com/science/article/pii/S0925400522017439#sec0055
                                   ion = c("neu"), # charge at pH 7.4
                                   VP=NA,
                                   molW=643.26,
                                   waterSol=NA) 

physchem_T2_Phase_I <- data.frame(logKolw = c(NA),
                                  fup = c(NA),
                                  fut = c(NA),
                                  logKow = c(-0.1), # Unknown, but expected to be lower than parent compound. Used a delta of 1 based on https://pmc.ncbi.nlm.nih.gov/articles/PMC10935128/pdf/main.pdf
                                  ion = c("neu"), # charge at pH 7.4
                                  VP=NA,
                                  molW=483.22,
                                  waterSol=NA) 

physchem_HT2_Phase_II <- data.frame(logKolw = c(NA),
                                    fup = c(NA),
                                    fut = c(NA),
                                    logKow = c(-1.6), # Unknown, but expected to be lower than parent compound. Used a delta of 2 based on https://www.sciencedirect.com/science/article/pii/S0925400522017439#sec0055
                                    ion = c("neu"), # charge at pH 7.4
                                    VP=NA,
                                    molW=601.25,
                                    waterSol=NA) 

physchem_HT2_Phase_I <- data.frame(logKolw = c(NA),
                                   fup = c(NA),
                                   fut = c(NA),
                                   logKow = c(-0.6), # Unknown, but expected to be lower than parent compound. Used a delta of 1 based on https://pmc.ncbi.nlm.nih.gov/articles/PMC10935128/pdf/main.pdf
                                   ion = c("neu"), # charge at pH 7.4
                                   VP=NA,
                                   molW=441.21,
                                   waterSol=NA) 

# Exposure regimen ----
dose <- list()
for(i in 1:n_samples){ # length(IDs)){
  Total_exposure_select <- Total_exposure # %>% dplyr::filter(ID == IDs[i]) # use filter in case of >1 IDs
  T2_exposure_select <- Total_exposure_select %>% dplyr::filter(toxin == "T2")
  HT2_exposure_select <- Total_exposure_select %>% dplyr::filter(toxin == "HT2")
  
  T2_exposure_PBK <- T2_exposure_select %>% dplyr::select(Time_PBK, mean) %>%
    dplyr::rename("Time" = "Time_PBK") %>%
    dplyr::rename("exposure" = "mean") %>%
    dplyr::mutate("Time_adj" = NA)
  HT2_exposure_PBK <- HT2_exposure_select %>% dplyr::select(Time_PBK, mean) %>%
    dplyr::rename("Time" = "Time_PBK") %>%
    dplyr::rename("exposure" = "mean") %>%
    dplyr::mutate("Time_adj" = NA)
  
  # Adjust the time points of the exposure slightly so that they match with the time points for the simulation
  # otherwise, the Aoral_total is off
  # e.g. and exposure time point of 1.083 h will be adjusted to 1.1 h
  for(r in 1:nrow(T2_exposure_PBK)){
    T2_exposure_PBK$Time_adj[r] <- times[which.min(abs(times - T2_exposure_PBK$Time[r]))]  
  }
  
  for(r in 1:nrow(HT2_exposure_PBK)){
    HT2_exposure_PBK$Time_adj[r] <- times[which.min(abs(times - HT2_exposure_PBK$Time[r]))]  
  }
  
  T2_exposure_PBK <- T2_exposure_PBK %>% dplyr::select(Time_adj, exposure) %>%
    dplyr::rename("Time" = "Time_adj")
  
  HT2_exposure_PBK <- HT2_exposure_PBK %>% dplyr::select(Time_adj, exposure) %>%
    dplyr::rename("Time" = "Time_adj")
  
  regimen <- data.frame(Time = times, exposure = rep(0, length(times)))
  regimen_T2 <- rbind(regimen,T2_exposure_PBK)
  regimen_T2 <- dplyr::distinct(regimen_T2, Time)
  regimen_T2 <- left_join(regimen_T2,T2_exposure_PBK, by="Time") %>%
    mutate_if(is.numeric,coalesce,0)
  regimen_T2 <- dplyr::arrange(regimen_T2, Time)
  exposure_T2 <- aggregate(exposure~Time,regimen_T2,sum)
  exposure_T2$dose <- exposure_T2$exposure*n_step # to correct for the time steps
  exposure_T2$dose[1] <- exposure_T2$dose[1]*2 # to have a correct dose for the first time point
  
  regimen_HT2 <- rbind(regimen,HT2_exposure_PBK)
  regimen_HT2 <- dplyr::distinct(regimen_HT2, Time)
  regimen_HT2 <- left_join(regimen_HT2,HT2_exposure_PBK, by="Time") %>%
    mutate_if(is.numeric,coalesce,0)
  regimen_HT2 <- dplyr::arrange(regimen_HT2, Time)
  exposure_HT2 <- aggregate(exposure~Time,regimen_HT2,sum)
  exposure_HT2$dose <- exposure_HT2$exposure*n_step # to correct for the time steps
  exposure_HT2$dose[1] <- exposure_HT2$dose[1]*2 # to have a correct dose for the first time point
  
  dosing_T2 <- approxfun(x = exposure_T2$Time, y = exposure_T2$dose, rule = 2)
  dosing_T2 <- list(dosing = dosing_T2)
  dosing_HT2 <- approxfun(x = exposure_HT2$Time, y = exposure_HT2$dose, rule = 2)
  dosing_HT2 <- list(dosing = dosing_HT2)
  
  dose[[i]] <- IDs[i]
  dose[[i]][[1]] <- list(exposure_T2$Time)
  dose[[i]][[2]] <- dosing_T2
  dose[[i]][[3]] <- dosing_HT2
}

# Set initial values ----
A.ini_T2_HT2 <- c(Aalv_T2 = 0, Alun_T2 = 0, Ass_T2 = 0, AskE_T2 = 0, AskU_T2 = 0, 
                  Aoral_total_T2 = 0, AstoL_T2 = 0, AintL_T2 = 0, Aint_metab_T2 = 0,  Asto_T2 = 0, Aint_T2 = 0, Afec_T2 = 0, 
                  Ahrt_T2 = 0, Abrn_T2 = 0, Athy_T2 = 0, Aadp_T2 = 0, Amus_T2 = 0, Abon_T2 = 0, Agon_T2 = 0, Arem_T2 = 0,
                  Apan_T2 = 0, Aspl_T2 = 0, Aliv_T2 = 0, AlivCL_T2 = 0, Akid_T2 = 0, AkidCL_T2 = 0,
                  Aa_T2 = 0, Av_T2 = 0, CumulativeIntake_T2 = 0,
                  
                  Aalv_HT2 = 0, Alun_HT2 = 0, Ass_HT2 = 0, AskE_HT2 = 0, AskU_HT2 = 0, 
                  Aoral_total_HT2 = 0, AstoL_HT2 = 0, AintL_HT2 = 0, Aint_metab_HT2 = 0, Asto_HT2 = 0, Aint_HT2 = 0, Afec_HT2 = 0, 
                  Ahrt_HT2 = 0, Abrn_HT2 = 0, Athy_HT2 = 0, Aadp_HT2 = 0, Amus_HT2 = 0, Abon_HT2 = 0, Agon_HT2 = 0, Arem_HT2 = 0,
                  Apan_HT2 = 0, Aspl_HT2 = 0, Aliv_HT2 = 0, AlivCL_HT2 = 0, Akid_HT2 = 0, AkidCL_HT2 = 0,
                  Aa_HT2 = 0, Av_HT2 = 0, CumulativeIntake_HT2 = 0,
                  
                  Aalv_T2_Phase_I = 0, Alun_T2_Phase_I = 0, Ass_T2_Phase_I = 0, AskE_T2_Phase_I = 0, AskU_T2_Phase_I = 0, 
                  Aoral_total_T2_Phase_I = 0, AstoL_T2_Phase_I = 0, AintL_T2_Phase_I = 0, Aint_metab_T2_Phase_I = 0, Asto_T2_Phase_I = 0, Aint_T2_Phase_I = 0, Afec_T2_Phase_I = 0, 
                  Ahrt_T2_Phase_I = 0, Abrn_T2_Phase_I = 0, Athy_T2_Phase_I = 0, Aadp_T2_Phase_I = 0, Amus_T2_Phase_I = 0, Abon_T2_Phase_I = 0, Agon_T2_Phase_I = 0, Arem_T2_Phase_I = 0,
                  Apan_T2_Phase_I = 0, Aspl_T2_Phase_I = 0, Aliv_T2_Phase_I = 0, AlivCL_T2_Phase_I = 0, Akid_T2_Phase_I = 0, AkidCL_T2_Phase_I = 0,
                  Aa_T2_Phase_I = 0, Av_T2_Phase_I = 0, CumulativeIntake_T2_Phase_I = 0,
                  
                  Aalv_HT2_Phase_I = 0, Alun_HT2_Phase_I = 0, Ass_HT2_Phase_I = 0, AskE_HT2_Phase_I = 0, AskU_HT2_Phase_I = 0, 
                  Aoral_total_HT2_Phase_I = 0, AstoL_HT2_Phase_I = 0, AintL_HT2_Phase_I = 0, Aint_metab_HT2_Phase_I = 0, Asto_HT2_Phase_I = 0, Aint_HT2_Phase_I = 0, Afec_HT2_Phase_I = 0, 
                  Ahrt_HT2_Phase_I = 0, Abrn_HT2_Phase_I = 0, Athy_HT2_Phase_I = 0, Aadp_HT2_Phase_I = 0, Amus_HT2_Phase_I = 0, Abon_HT2_Phase_I = 0, Agon_HT2_Phase_I = 0, Arem_HT2_Phase_I = 0,
                  Apan_HT2_Phase_I = 0, Aspl_HT2_Phase_I = 0, Aliv_HT2_Phase_I = 0, AlivCL_HT2_Phase_I = 0, Akid_HT2_Phase_I = 0, AkidCL_HT2_Phase_I = 0,
                  Aa_HT2_Phase_I = 0, Av_HT2_Phase_I = 0, CumulativeIntake_HT2_Phase_I = 0,
                  
                  Aalv_T2_Phase_II = 0, Alun_T2_Phase_II = 0, Ass_T2_Phase_II = 0, AskE_T2_Phase_II = 0, AskU_T2_Phase_II = 0, 
                  Aoral_total_T2_Phase_II = 0, AstoL_T2_Phase_II = 0, AintL_T2_Phase_II = 0, Aint_metab_T2_Phase_II = 0, Asto_T2_Phase_II = 0, Aint_T2_Phase_II = 0, Afec_T2_Phase_II = 0, 
                  Ahrt_T2_Phase_II = 0, Abrn_T2_Phase_II = 0, Athy_T2_Phase_II = 0, Aadp_T2_Phase_II = 0, Amus_T2_Phase_II = 0, Abon_T2_Phase_II = 0, Agon_T2_Phase_II = 0, Arem_T2_Phase_II = 0,
                  Apan_T2_Phase_II = 0, Aspl_T2_Phase_II = 0, Aliv_T2_Phase_II = 0, AlivCL_T2_Phase_II = 0, Akid_T2_Phase_II = 0, AkidCL_T2_Phase_II = 0,
                  Aa_T2_Phase_II = 0, Av_T2_Phase_II = 0, CumulativeIntake_T2_Phase_II = 0,
                  
                  Aalv_HT2_Phase_II = 0, Alun_HT2_Phase_II = 0, Ass_HT2_Phase_II = 0, AskE_HT2_Phase_II = 0, AskU_HT2_Phase_II = 0, 
                  Aoral_total_HT2_Phase_II = 0, AstoL_HT2_Phase_II = 0, AintL_HT2_Phase_II = 0, Aint_metab_HT2_Phase_II = 0, Asto_HT2_Phase_II = 0, Aint_HT2_Phase_II = 0, Afec_HT2_Phase_II = 0, 
                  Ahrt_HT2_Phase_II = 0, Abrn_HT2_Phase_II = 0, Athy_HT2_Phase_II = 0, Aadp_HT2_Phase_II = 0, Amus_HT2_Phase_II = 0, Abon_HT2_Phase_II = 0, Agon_HT2_Phase_II = 0, Arem_HT2_Phase_II = 0,
                  Apan_HT2_Phase_II = 0, Aspl_HT2_Phase_II = 0, Aliv_HT2_Phase_II = 0, AlivCL_HT2_Phase_II = 0, Akid_HT2_Phase_II = 0, AkidCL_HT2_Phase_II = 0,
                  Aa_HT2_Phase_II = 0, Av_HT2_Phase_II = 0, CumulativeIntake_HT2_Phase_II = 0)

# Define activity pattern ---- 
# In this case, individual is at rest
activity_pattern <- data.frame( level   = c('light', 'moderate', 'heavy'), # Activity level.  
                                n_day   = c(0, 0,  0), # number of activities per day
                                start   = c(6, 12, 18), # Start time (h) of the exercise level on a day
                                dur     = c(6, 6,  6), # duration of the activity
                                int     = c(24, 24, 24), # if multiple events take place, interval between the start of the activities. (e.g. if an activity start at 8h, with an interval of 4, then the second activity starts at 12h)
                                t_stop  = c(1*24,1*24, 1*24)) # time (d) when activity a particular activity level should stop.

activity_pattern <- f.activityreg(activity_pattern, sim_time, n_step)

# Define kinetic parameters ----
kinetic_parameters <- data.frame(
  fss = rep(0.01,n_samples),
  Vss = rep(1,n_samples), 
  kse    = rlnorm(n_samples,mean=log(4),sd=0.2), # stomach emptying rate constant (h^-1)
  kfec   = rlnorm(n_samples,mean=log(0.3),sd=0.1), # fecal excretion rate constant (h^-1)
  
  ka_sto_T2   = rnorm(n_samples,mean=ka_T2,sd=(0.2*ka_T2)), # absorption rate constant from stomach lumen (h^-1)
  ka_int_T2   = rnorm(n_samples,mean=ka_T2,sd=(0.2*ka_T2)), # absorption rate constant from intestinal lumen (h^-1)
  km_sto_T2_HT2 = rnorm(n_samples,mean=km_sto_T2_HT2,sd=(0.2*km_sto_T2_HT2)),
  km_sto_T2_T2_Phase_II = rnorm(n_samples,mean=km_sto_T2_T2_Phase_II,sd=(0.2*km_sto_T2_T2_Phase_II)),
  km_sto_T2_T2_Phase_I = rnorm(n_samples,mean=km_sto_T2_T2_Phase_I,sd=(0.2*km_sto_T2_T2_Phase_I)),
  km_int_T2_HT2 = rnorm(n_samples,mean=km_int_T2_HT2,sd=(0.2*km_int_T2_HT2)),
  km_int_T2_T2_Phase_II = rnorm(n_samples,mean=km_int_T2_T2_Phase_II,sd=(0.2*km_int_T2_T2_Phase_II)),
  km_int_T2_T2_Phase_I = rnorm(n_samples,mean=km_int_T2_T2_Phase_I,sd=(0.2*km_int_T2_T2_Phase_I)),
  CLliv_T2_HT2 = rnorm(n_samples,mean=CLliv_T2_HT2,sd=(0.2*CLliv_T2_HT2)),
  CLliv_T2_T2_Phase_II = rnorm(n_samples,mean=CLliv_T2_T2_Phase_II,sd=(0.2*CLliv_T2_T2_Phase_II)),
  CLliv_T2_T2_Phase_I = rnorm(n_samples,mean=CLliv_T2_T2_Phase_I,sd=(0.2*CLliv_T2_T2_Phase_I)),
  CLkid_T2 = rnorm(n_samples,mean=CLkid_T2,sd=(0.2*CLkid_T2)),
  PAalv_T2 = rep(0,n_samples),
  Pluna_T2  = rep(1E4,n_samples),
  Ks_T2 = rep(1E-4,n_samples),
  
  ka_sto_HT2 = rnorm(n_samples,mean=ka_HT2,sd=(0.2*ka_HT2)), # absorption rate constant from stomach lumen (h^-1)
  ka_int_HT2 = rnorm(n_samples,mean=ka_HT2,sd=(0.2*ka_HT2)), # absorption rate constant from intestinal lumen (h^-1)
  km_sto_HT2_HT2_Phase_II = rnorm(n_samples,mean=km_sto_HT2_HT2_Phase_II,sd=(0.2*km_sto_HT2_HT2_Phase_II)),
  km_sto_HT2_HT2_Phase_I = rnorm(n_samples,mean=km_sto_HT2_HT2_Phase_I,sd=(0.2*km_sto_HT2_HT2_Phase_I)),
  km_int_HT2_HT2_Phase_II = rnorm(n_samples,mean=km_int_HT2_HT2_Phase_II,sd=(0.2*km_int_HT2_HT2_Phase_II)),
  km_int_HT2_HT2_Phase_I = rnorm(n_samples,mean=km_int_HT2_HT2_Phase_I,sd=(0.2*km_int_HT2_HT2_Phase_I)),
  CLliv_HT2_HT2_Phase_II = rnorm(n_samples,mean=CLliv_HT2_HT2_Phase_II,sd=(0.2*CLliv_HT2_HT2_Phase_II)),
  CLliv_HT2_HT2_Phase_I = rnorm(n_samples,mean=CLliv_HT2_HT2_Phase_I,sd=(0.2*CLliv_HT2_HT2_Phase_I)),
  CLkid_HT2 = rnorm(n_samples,mean=CLkid_HT2,sd=(0.2*CLkid_HT2)),
  PAalv_HT2 = rep(0,n_samples),
  Pluna_HT2  = rep(1E4,n_samples),
  Ks_HT2 = rep(1E-4,n_samples),
  
  ka_sto_T2_Phase_I = rnorm(n_samples,mean=ka_T2,sd=(0.2*ka_T2)), # absorption rate constant from stomach lumen (h^-1)
  ka_int_T2_Phase_I = rnorm(n_samples,mean=ka_T2,sd=(0.2*ka_T2)), # absorption rate constant from intestinal lumen (h^-1)
  CLkid_T2_Phase_I = rnorm(n_samples,mean=CLkid_T2_Phase_I,sd=(0.2*CLkid_T2_Phase_I)),
  PAalv_T2_Phase_I = rep(0,n_samples),
  Pluna_T2_Phase_I  = rep(1E4,n_samples),
  Ks_T2_Phase_I = rep(1E-4,n_samples),
  
  ka_sto_HT2_Phase_I = rnorm(n_samples,mean=ka_HT2,sd=(0.2*ka_HT2)), # absorption rate constant from stomach lumen (h^-1)
  ka_int_HT2_Phase_I = rnorm(n_samples,mean=ka_HT2,sd=(0.2*ka_HT2)), # absorption rate constant from intestinal lumen (h^-1)
  CLkid_HT2_Phase_I = rnorm(n_samples,mean=CLkid_HT2_Phase_I,sd=(0.2*CLkid_HT2_Phase_I)),
  PAalv_HT2_Phase_I = rep(0,n_samples),
  Pluna_HT2_Phase_I  = rep(1E4,n_samples),
  Ks_HT2_Phase_I = rep(1E-4,n_samples),
  
  ka_sto_T2_Phase_II = rnorm(n_samples,mean=ka_T2,sd=(0.2*ka_T2)), # absorption rate constant from stomach lumen (h^-1)
  ka_int_T2_Phase_II = rnorm(n_samples,mean=ka_T2,sd=(0.2*ka_T2)), # absorption rate constant from intestinal lumen (h^-1)
  km_sto_T2_Phase_II_HT2_Phase_II = rnorm(n_samples,mean=km_sto_T2_Phase_II_HT2_Phase_II,sd=(0.2*km_sto_T2_Phase_II_HT2_Phase_II)),
  km_int_T2_Phase_II_HT2_Phase_II = rnorm(n_samples,mean=km_int_T2_Phase_II_HT2_Phase_II,sd=(0.2*km_int_T2_Phase_II_HT2_Phase_II)),
  CLliv_T2_Phase_II_HT2_Phase_II = rnorm(n_samples,mean=CLliv_T2_Phase_II_HT2_Phase_II,sd=(0.2*CLliv_T2_Phase_II_HT2_Phase_II)),
  CLkid_T2_Phase_II = rnorm(n_samples,mean=CLkid_T2_Phase_II,sd=(0.2*CLkid_T2_Phase_II)),
  PAalv_T2_Phase_II = rep(0,n_samples),
  Pluna_T2_Phase_II  = rep(1E4,n_samples),
  Ks_T2_Phase_II = rep(1E-4,n_samples),
  
  ka_sto_HT2_Phase_II = rnorm(n_samples,mean=ka_HT2,sd=(0.2*ka_HT2)), # absorption rate constant from stomach lumen (h^-1)
  ka_int_HT2_Phase_II = rnorm(n_samples,mean=ka_HT2,sd=(0.2*ka_HT2)), # absorption rate constant from intestinal lumen (h^-1)
  CLkid_HT2_Phase_II = rnorm(n_samples,mean=CLkid_HT2_Phase_II,sd=(0.2*CLkid_HT2_Phase_II)),
  PAalv_HT2_Phase_II = rep(0,n_samples),
  Pluna_HT2_Phase_II  = rep(1E4,n_samples),
  Ks_HT2_Phase_II = rep(1E-4,n_samples)
)

# Run model ----
out <- list()
physB_overview <- data.frame()
for (n in 1:n_samples) {
  print(paste("Simulating subject ",n," of ",n_samples,sep=""))
  physiology_df_select <- physiology_df[n,]
  physB <- f.physB(physiology_df_select)
  #  Find and replace negative adipose tissue compartment volumes
  while (physB$Vadp < 0) {
    physB <- f.physB(physiology_df_select)
  }
  physB_overview <- rbind(physB_overview,physB)
  
  Ptisb_T2 <- get_tissue_blood_partition(physchem_T2$VP, physchem_T2$fup, physchem_T2$fut, physchem_T2$logKow, physchem_T2$ion, physchem_T2$VP, physchem_T2$molW, physchem_T2$waterSol)
  Ptisb_HT2 <- get_tissue_blood_partition(physchem_HT2$VP, physchem_HT2$fup, physchem_HT2$fut, physchem_HT2$logKow, physchem_HT2$ion, physchem_HT2$VP, physchem_HT2$molW, physchem_HT2$waterSol)
  Ptisb_T2_Phase_II <- get_tissue_blood_partition(physchem_T2_Phase_II$VP, physchem_T2_Phase_II$fup, physchem_T2_Phase_II$fut, physchem_T2_Phase_II$logKow, physchem_T2_Phase_II$ion, physchem_T2_Phase_II$VP, physchem_T2_Phase_II$molW, physchem_T2_Phase_II$waterSol)
  Ptisb_HT2_Phase_II <- get_tissue_blood_partition(physchem_HT2_Phase_II$VP, physchem_HT2_Phase_II$fup, physchem_HT2_Phase_II$fut, physchem_HT2_Phase_II$logKow, physchem_HT2_Phase_II$ion, physchem_HT2_Phase_II$VP, physchem_HT2_Phase_II$molW, physchem_HT2_Phase_II$waterSol)
  Ptisb_T2_Phase_I <- get_tissue_blood_partition(physchem_T2_Phase_I$VP, physchem_T2_Phase_I$fup, physchem_T2_Phase_I$fut, physchem_T2_Phase_I$logKow, physchem_T2_Phase_I$ion, physchem_T2_Phase_I$VP, physchem_T2_Phase_I$molW, physchem_T2_Phase_I$waterSol)
  Ptisb_HT2_Phase_I <- get_tissue_blood_partition(physchem_HT2_Phase_I$VP, physchem_HT2_Phase_I$fup, physchem_HT2_Phase_I$fut, physchem_HT2_Phase_I$logKow, physchem_HT2_Phase_I$ion, physchem_HT2_Phase_I$VP, physchem_HT2_Phase_I$molW, physchem_HT2_Phase_I$waterSol)
  
  parms <- f.parms.pbk(physB, Ptisb_T2, Ptisb_HT2, Ptisb_T2_Phase_II, Ptisb_HT2_Phase_II, Ptisb_T2_Phase_I, Ptisb_HT2_Phase_I)
  
  parms["Ks_T2"] <- kinetic_parameters$Ks_T2[n]
  parms["fss"] <- kinetic_parameters$fss[n]
  parms["Vss"] <- kinetic_parameters$Vss[n]
  parms["kse"] <- kinetic_parameters$kse[n]
  parms["kfec"] <- kinetic_parameters$kfec[n]
  
  parms["ka_sto_T2"] <- kinetic_parameters$ka_sto_T2[n]
  parms["ka_int_T2"] <- kinetic_parameters$ka_int_T2[n]
  parms["km_sto_T2_HT2"] <- kinetic_parameters$km_sto_T2_HT2[n]
  parms["km_sto_T2_T2_Phase_II"] <- kinetic_parameters$km_sto_T2_T2_Phase_II[n]
  parms["km_sto_T2_T2_Phase_I"] <- kinetic_parameters$km_sto_T2_T2_Phase_I[n]
  parms["km_int_T2_HT2"] <- kinetic_parameters$km_int_T2_HT2[n]
  parms["km_int_T2_T2_Phase_II"] <- kinetic_parameters$km_int_T2_T2_Phase_II[n]
  parms["km_int_T2_T2_Phase_I"] <- kinetic_parameters$km_int_T2_T2_Phase_I[n]
  parms["CLliv_T2_HT2"] <- kinetic_parameters$CLliv_T2_HT2[n]
  parms["CLliv_T2_T2_Phase_II"] <- kinetic_parameters$CLliv_T2_T2_Phase_II[n]
  parms["CLliv_T2_T2_Phase_I"] <- kinetic_parameters$CLliv_T2_T2_Phase_I[n]
  parms["CLkid_T2"] <- kinetic_parameters$CLkid_T2[n]
  
  parms["ka_sto_HT2"] <- kinetic_parameters$ka_sto_HT2[n]
  parms["ka_int_HT2"] <- kinetic_parameters$ka_int_HT2[n]
  parms["km_sto_HT2_HT2_Phase_II"] <- kinetic_parameters$km_sto_HT2_HT2_Phase_II[n]
  parms["km_sto_HT2_HT2_Phase_I"] <- kinetic_parameters$km_sto_HT2_HT2_Phase_I[n]
  parms["km_int_HT2_HT2_Phase_II"] <- kinetic_parameters$km_int_HT2_HT2_Phase_II[n]
  parms["km_int_HT2_HT2_Phase_I"] <- kinetic_parameters$km_int_HT2_HT2_Phase_I[n]
  parms["CLliv_HT2_HT2_Phase_II"] <- kinetic_parameters$CLliv_HT2_HT2_Phase_II[n]
  parms["CLliv_HT2_HT2_Phase_I"] <- kinetic_parameters$CLliv_HT2_HT2_Phase_I[n]
  parms["CLkid_HT2"] <- kinetic_parameters$CLkid_HT2[n]
  
  parms["ka_sto_T2_Phase_I"] <- kinetic_parameters$ka_sto_T2_Phase_I[n]
  parms["ka_int_T2_Phase_I"] <- kinetic_parameters$ka_int_T2_Phase_I[n]
  parms["km_sto_T2_Phase_I"] <- kinetic_parameters$km_sto_T2_Phase_I[n]
  parms["km_int_T2_Phase_I"] <- kinetic_parameters$km_int_T2_Phase_I[n]
  parms["CLliv_T2_Phase_I"] <- kinetic_parameters$CLliv_T2_Phase_I[n]
  parms["CLkid_T2_Phase_I"] <- kinetic_parameters$CLkid_T2_Phase_I[n]
  
  parms["ka_sto_HT2_Phase_I"] <- kinetic_parameters$ka_sto_HT2_Phase_I[n]
  parms["ka_int_HT2_Phase_I"] <- kinetic_parameters$ka_int_HT2_Phase_I[n]
  parms["CLkid_HT2_Phase_I"] <- kinetic_parameters$CLkid_HT2_Phase_I[n]
  
  parms["ka_sto_T2_Phase_II"] <- kinetic_parameters$ka_sto_T2_Phase_II[n]
  parms["ka_int_T2_Phase_II"] <- kinetic_parameters$ka_int_T2_Phase_II[n]
  parms["km_sto_T2_Phase_II_HT2_Phase_II"] <- kinetic_parameters$km_sto_T2_Phase_II_HT2_Phase_II[n]
  parms["km_int_T2_Phase_II_HT2_Phase_II"] <- kinetic_parameters$km_int_T2_Phase_II_HT2_Phase_II[n]
  parms["CLliv_T2_Phase_II_HT2_Phase_II"] <- kinetic_parameters$CLliv_T2_Phase_II_HT2_Phase_II[n]
  parms["CLkid_T2_Phase_II"] <- kinetic_parameters$CLkid_T2_Phase_II[n]
  
  parms["ka_sto_HT2_Phase_II"] <- kinetic_parameters$ka_sto_HT2_Phase_II[n]
  parms["ka_int_HT2_Phase_II"] <- kinetic_parameters$ka_int_HT2_Phase_II[n]
  parms["CLkid_HT2_Phase_II"] <- kinetic_parameters$CLkid_HT2_Phase_II[n]
  
  if (include_lung==TRUE) {
    parms["PAalv_T2"] <- kinetic_parameters$PAalv_T2[n]
    parms["Pluna_T2"] <- kinetic_parameters$Pluna_T2[n]
    parms["PAalv_HT2"] <- kinetic_parameters$PAalv_HT2[n]
    parms["Pluna_HT2"] <- kinetic_parameters$Pluna_HT2[n]
    parms["PAalv_T2_Phase_I"] <- kinetic_parameters$PAalv_T2_Phase_I[n]
    parms["Pluna_T2_Phase_I"] <- kinetic_parameters$Pluna_T2_Phase_I[n]
    parms["PAalv_HT2_Phase_I"] <- kinetic_parameters$PAalv_HT2_Phase_I[n]
    parms["Pluna_HT2_Phase_I"] <- kinetic_parameters$Pluna_HT2_Phase_I[n]
    parms["PAalv_T2_Phase_II"] <- kinetic_parameters$PAalv_T2_Phase_II[n]
    parms["Pluna_T2_Phase_II"] <- kinetic_parameters$Pluna_T2_Phase_II[n]
    parms["PAalv_HT2_Phase_II"] <- kinetic_parameters$PAalv_HT2_Phase_II[n]
    parms["Pluna_HT2_Phase_II"] <- kinetic_parameters$Pluna_HT2_Phase_II[n]
    
  } else {
    parms["Pab_T2"] <- kinetic_parameters$Pab_T2[n]
    parms["Pab_HT2"] <- kinetic_parameters$Pab_HT2[n]
    parms["Pab_T2_Phase_I"] <- kinetic_parameters$Pab_T2_Phase_I[n]
    parms["Pab_HT2_Phase_I"] <- kinetic_parameters$Pab_HT2_Phase_I[n]
    parms["Pab_T2_Phase_II"] <- kinetic_parameters$Pab_T2_Phase_II[n]
    parms["Pab_HT2_Phase_II"] <- kinetic_parameters$Pab_HT2_Phase_II[n]
  }
  # Correct parms for activity levels
  parms = f.activity(parms, activity_pattern, sex=physiology_df_select$Sex, approx_fun=FALSE)
  
  times_PBK <- dose[[n]][[1]][[1]]
  R.se   <- rep(parms$kse, length(times_PBK))
  
  # Dosing to compartments
  X.oral_T2 <- dose[[n]][[2]]$dosing #
  X.oral_HT2 <- dose[[n]][[3]]$dosing #
  
  out[[n]] <- eval(as.data.frame(lsoda(A.ini_T2_HT2, times_PBK, f.pbk, parms = parms,
                                       R.se = R.se, tstep = 1/n_step,
                                       f.select = f.select, include_lung = include_lung,
                                       events=list(data=events))))
}

# Save output ----
save(out, file=paste("output_T2_HT2_Phase_I_II_",Sys.Date(),".Rdata",sep=""))

# Check output ----
for(n in 1:n_samples){
  output_n <- out[[n]]
  print(max(output_n$Aurine_HT2))
}

for(n in 1:n_samples){
  output_n <- out[[n]]
  print(max(output_n$Ca_HT2))
}

# Figures ----
myColors <- paletteer_c("grDevices::rainbow", n_samples)

## Individual plot ---- 
for(n in 1:n_samples){
  Urinedata_select <- Urinedata # %>% dplyr::filter(ID == IDs[n]) # use filter in case of >1 IDs
  
  output_n <- out[[n]]
  Figure_Aurine <- ggplot()
  Figure_Aurine <- Figure_Aurine +
    geom_line(data=output_n, aes(x=time, y=CumulativeIntake_T2, color="001_T2"), lty='dashed') +
    geom_line(data=output_n, aes(x=time, y=CumulativeIntake_HT2, color="002_HT2"), lty='dashed') +
    geom_line(data=output_n, aes(x=time, y=AkidCL_T2, color="01_T2")) +
    geom_line(data=output_n, aes(x=time, y=AkidCL_T2_Phase_I, color="02_T2_Phase_I")) +
    geom_line(data=output_n, aes(x=time, y=AkidCL_T2_Phase_II, color="03_T2_Phase_II")) +
    geom_line(data=output_n, aes(x=time, y=Aurine_T2, color="04_T2_Urine")) +
    geom_line(data=output_n, aes(x=time, y=AkidCL_HT2, color="05_HT2")) +
    geom_line(data=output_n, aes(x=time, y=AkidCL_HT2_Phase_I, color="06_HT2_Phase_I")) +
    geom_line(data=output_n, aes(x=time, y=AkidCL_HT2_Phase_II, color="07_HT2_Phase_II")) +
    geom_line(data=output_n, aes(x=time, y=Aurine_HT2, color="08_HT2_Urine")) +
    geom_point(data=Urinedata_select, aes(x=Time_PBK, y=Cumulative_amount_HT2, color="09_HT2_Phase_II")) +
    
    scale_colour_manual(name='',
                        values=c('001_T2'='red',
                                 '002_HT2'='blue',
                                 '01_T2'='red',
                                 '02_T2_Phase_I'='magenta',
                                 '03_T2_Phase_II'='palevioletred',
                                 '04_T2_Urine'='firebrick4',
                                 '05_HT2'='blue',
                                 '06_HT2_Phase_I'='cyan',
                                 '07_HT2_Phase_II'='dodgerblue',
                                 '08_HT2_Urine'='darkblue',
                                 '09_HT2_Phase_II'='darkblue'),
                        labels=c('001_T2'='Total T-2 intake',
                                 '002_HT2'='Total HT-2 intake',
                                 '01_T2'='Unchanged T-2 in urine',
                                 '02_T2_Phase_I'='Phase I metabolites of T-2 in urine',
                                 '03_T2_Phase_II'='Phase II metabolites of T-2 in urine',
                                 '04_T2_Urine'='Total unchanged T-2 and phase II metabolites of T-2 in urine',
                                 '05_HT2'='Unchanged HT-2 in urine',
                                 '06_HT2_Phase_I'='Phase I metabolites of HT-2 in urine',
                                 '07_HT2_Phase_II'='Phase II metabolites of HT-2 in urine',
                                 '08_HT2_Urine'='Total unchanged HT-2 and phase II metabolites of HT-2 in urine',
                                 '09_HT2_Phase_II'='Total unchanged HT-2 and phase II metabolites of HT-2 in urine')) +
    
    theme_bw() +
    labs(title=paste("Cumulative amounts of T-2/HT-2 and their metabolites - Subject ",n,sep=""),x="\nTime (h)", y="Amount (ng)\n") + # \u03BC = micro symbol
    scale_x_continuous(breaks=c(0,4,8,12,16,20,24), limits=c(0,24)) +
    theme(plot.title = element_text(hjust = 0.5))
  print(Figure_Aurine)
}

## Plot with CI ----
output_sim <- data.frame()
for(i in 1:n_samples){
  output_n <- out[[i]]
  output_n$Subject <- i
  output_sim <- rbind(output_sim,output_n)
}

output_sim_wide <- reshape(output_sim, idvar = "time", timevar = "Subject", direction = "wide")
output_sim_AkidCL_T2 <- output_sim_wide %>%
  dplyr::select(time,starts_with("AkidCL_T2."))

output_sim_AkidCL_T2_Phase_I <- output_sim_wide %>%
  dplyr::select(time,starts_with("AkidCL_T2_Phase_I."))

output_sim_AkidCL_T2_Phase_II <- output_sim_wide %>%
  dplyr::select(time,starts_with("AkidCL_T2_Phase_II."))

output_sim_AkidCL_HT2 <- output_sim_wide %>%
  dplyr::select(time,starts_with("AkidCL_HT2."))

output_sim_AkidCL_HT2_Phase_I <- output_sim_wide %>%
  dplyr::select(time,starts_with("AkidCL_HT2_Phase_I."))

output_sim_AkidCL_HT2_Phase_II <- output_sim_wide %>%
  dplyr::select(time,starts_with("AkidCL_HT2_Phase_II."))

output_sim_AkidCL_T2_range <- cbind(
  time = times,
  as.data.frame(t(apply(output_sim_AkidCL_T2[-1], 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T))))))

output_sim_AkidCL_T2_Phase_I_range <- cbind(
  time = times,
  as.data.frame(t(apply(output_sim_AkidCL_T2_Phase_I[-1], 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T))))))

output_sim_AkidCL_T2_Phase_II_range <- cbind(
  time = times,
  as.data.frame(t(apply(output_sim_AkidCL_T2_Phase_II[-1], 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T))))))

output_sim_AkidCL_HT2_range <- cbind(
  time = times,
  as.data.frame(t(apply(output_sim_AkidCL_HT2[-1], 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T))))))

output_sim_AkidCL_HT2_Phase_I_range <- cbind(
  time = times,
  as.data.frame(t(apply(output_sim_AkidCL_HT2_Phase_I[-1], 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T))))))

output_sim_AkidCL_HT2_Phase_II_range <- cbind(
  time = times,
  as.data.frame(t(apply(output_sim_AkidCL_HT2_Phase_II[-1], 1, function(y_est) c(
    median_est         = median(y_est,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T))))))

CI.plot <- ggplot() +
  geom_ribbon(data = output_sim_AkidCL_T2_range, aes(x = time, ymin = ci_lower_est, ymax = ci_upper_est), fill="red", alpha=0.3) +
  geom_line(data= output_sim_AkidCL_T2_range, aes(x = time, y = median_est, colour = "01_T2")) +
  geom_ribbon(data = output_sim_AkidCL_T2_Phase_I_range, aes(x = time, ymin = ci_lower_est, ymax = ci_upper_est), fill="magenta", alpha=0.3) +
  geom_line(data= output_sim_AkidCL_T2_Phase_I_range, aes(x = time, y = median_est, colour = "02_T2_Phase_I")) +
  geom_ribbon(data = output_sim_AkidCL_T2_Phase_II_range, aes(x = time, ymin = ci_lower_est, ymax = ci_upper_est), fill="palevioletred", alpha=0.3) +
  geom_line(data= output_sim_AkidCL_T2_Phase_II_range, aes(x = time, y = median_est, colour = "03_T2_Phase_II")) +
  geom_ribbon(data = output_sim_AkidCL_HT2_range, aes(x = time, ymin = ci_lower_est, ymax = ci_upper_est), fill="blue", alpha=0.3) +
  geom_line(data= output_sim_AkidCL_HT2_range, aes(x = time, y = median_est, colour = "04_HT2")) +
  geom_ribbon(data = output_sim_AkidCL_HT2_Phase_I_range, aes(x = time, ymin = ci_lower_est, ymax = ci_upper_est), fill="cyan", alpha=0.3) +
  geom_line(data= output_sim_AkidCL_HT2_Phase_I_range, aes(x = time, y = median_est, colour = "05_HT2_Phase_I")) +
  geom_ribbon(data = output_sim_AkidCL_HT2_Phase_II_range, aes(x = time, ymin = ci_lower_est, ymax = ci_upper_est), fill="dodgerblue", alpha=0.3) +
  geom_line(data= output_sim_AkidCL_HT2_Phase_II_range, aes(x = time, y = median_est, colour = "06_HT2_Phase_II")) +
  
  scale_colour_manual(name='Output',
                      values=c('01_T2'='red',
                               '02_T2_Phase_I'='magenta',
                               '03_T2_Phase_II'='palevioletred',
                               '04_HT2'='blue',
                               '05_HT2_Phase_I'='cyan',
                               '06_HT2_Phase_II'='dodgerblue'),
                      labels=c('01_T2'='Unchanged T-2 in urine',
                               '02_T2_Phase_I'='Phase I metabolites of T-2 in urine',
                               '03_T2_Phase_II'='Phase II metabolites of T-2 in urine',
                               '04_HT2'='Unchanged HT-2 in urine',
                               '05_HT2_Phase_I'='Phase I metabolites of HT-2 in urine',
                               '06_HT2_Phase_II'='Phase II metabolites of HT-2 in urine')) +
  
  ylab("Amount (ng)\n") + xlab("\nTime (h)") + theme_bw() +
  ggtitle("Cumulative amount of T-2/HT-2 and their metabolites in urine") +
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24), limits=c(0,24)) +
  theme(plot.title = element_text(hjust = 0.5))
CI.plot

# Calculating the fraction absorbed ----
Fabs <- data.frame("ID"=rep(NA,n_samples),
                   "T2_intake"=NA,
                   "HT2_intake"=NA,
                   "Total_intake"=NA,
                   "T2_urine"=NA,
                   "HT2_urine"=NA,
                   "Total_urine"=NA,
                   "HT2_gluc"=NA,
                   "Fabs"=NA)
for(n in 1:n_samples){
  Prediction <- out[[n]] 
  Time_last <- max(Prediction$time)
  Prediction <- out[[n]] %>%
    filter(time == Time_last) 
  
  Fabs$ID[n] <- IDs[n]
  Fabs$T2_intake[n] <- Prediction$CumulativeIntake_T2[1]
  Fabs$HT2_intake[n] <- Prediction$CumulativeIntake_HT2[1]
  Fabs$Total_intake[n] <- Fabs$T2_intake[n] + Fabs$HT2_intake[n]
  Fabs$T2_urine[n] <- Prediction$AkidCL_T2[1] + Prediction$AkidCL_T2_Phase_I[1] + Prediction$AkidCL_T2_Phase_II[1]
  Fabs$HT2_urine[n] <- Prediction$AkidCL_HT2[1] + Prediction$AkidCL_HT2_Phase_I[1] + Prediction$AkidCL_HT2_Phase_II[1]
  Fabs$Total_urine[n] <- Fabs$T2_urine[n] + Fabs$HT2_urine[n]
  Fabs$HT2_gluc[n] <- Prediction$AkidCL_HT2_Phase_II[1]
  Fabs$Fabs[n] <- Fabs$Total_urine[n]/Fabs$Total_intake[n]
}

mean(Fabs$Fabs)
median(Fabs$Fabs)
min(Fabs$Fabs)
max(Fabs$Fabs)

mean(Fabs$HT2_gluc)
median(Fabs$HT2_gluc)
min(Fabs$HT2_gluc)
max(Fabs$HT2_gluc)
