In this project we used the Norwegian EuroMix biomonitoring study for probabilistic external exposure estimates of the mycotoxin T-2 and HT-2 from the diet. The external exposure was be compared with measured T-2 and HT-2 metabolites in the urine through statistical modelling. 

The EuroMix study is a small (n=144) bionomitoring study performed as a part of the EU funded scientific project *European Test and Risk Assessment Strategies for Mixtures (EuroMix)*. Participants recorded and weighed their food consumption (weighed food record) for a 24h period for two non-consecutive days (2-3 weeks between) in a diary. All urine voids were collected and stored separately and combined into time-pools. You can read more about the Euromix study in Husoy et al. 2019 (https://www.sciencedirect.com/science/article/pii/S0160412019306944?via%3Dihub). 

T-2 toxin (T-2) and HT2 toxin (HT-2) are mycotoxins of the type A trichothecenes. These mycotoxines are produced by fungi of the type Fusarium genera. T-2 and HT-2 are mainly reported to contaminate cereals, especially oats. In Europe, oats and oat products are noted to be the most susceptible commodities to T-2 and HT-2 contamination. Oats grown in Scandinavia have reportedly high concentrations.Please be aware that part of the work on this repository have been provided by Rudolf Hoogenveen at RIVM. This is clearly indicated in the titles below. It's higly appreciated that Rudolf is willing to publish his work on my repository.

# First paper 

The results are published in "Exploring the Relationship between Daily Intake and Renal Excretion of the Mycotoxins T-2 and HT-2 Toxin in Humans", by Hannah P. McKeon, Rudolf Hoogenveen, Marta M. Sopel, Marloes A. A. Schepens, Marcel J. B. Mengelers, Annick D. van den Brand, Judith A. de Heer, Anne Lise Brantseater, Maria Kalyva, Trine Husoy. You can read the paper here (https://pubmed.ncbi.nlm.nih.gov/40306543/).

## Modelling of T2 and HT2 occurance data (Rudolf)

Distributional characteristics of T2/HT2-concentrations in food products are available from literature, see the Excel-file. Examples of these characteristics are mean values, median values, maximum values, fractions with measurement values (given LOQ value), etc. These characteristics were found for different food products from different countries from different publications etc. Since these characteristics are difficult to compare (e.g. how to fit a model on both maximum and mean value?), I developed a maximum likelihood method to estimate the model parameters.

### Files

THe document where the occurence modelling is described  
[Documentation & results.docx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Documentation%20%26%20results.docx)

The R code to run the occurance modeling  
[Occurence data modelling.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Code/Occurence%20data%20modelling.R)

The occurrence data from literature as input to the occurence modelling  
[Occurrence data.xlsx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Occurrence%20data.xlsx)  
[Overview occurrence T2HT2.xlsx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Overview%20occurrence%20T2HT2.xlsx)

The output from the occurance modelling should be like this file and this is also an input to the exposure modelling
[Results_occurance_T2HT2.csv](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Results_occurance_T2HT2.csv)


## Exposure to T-2 and HT-2 toxin from diet


The exposure assessment were performed per meal for 40 participants from the EuroMix study who had high consumption of oat and/or oat containing food items.The exposure estimate was based on the individual consumption data from the weighted 24-h dietary records and concentration data of T-2 and HT-2 in grain such as oat, rye and wheat. Concentration data was obtained from the litterature, and the concentrations found in whole grains were used to estimate T-2 or HT-2 concentrations in food items (e.g. bread, cereals and oat milk). These estimations were based on recipes or declaration of contents.

Individual data from the EuroMix study cannot be openly shared, and we therefore created a Dummy file for the food consumption. The resulting file will be a stochastic rearrangement of the data (including the ID numbers) in the original file, but all data will be maintained. The calculations with these Dummy data will therefore not any longer be done for individuals with high oat consumption. In addition the amounts will not match the right food category, and therefore the exposure from some food categories will not be reasonable. The Dummy file is only added for you to run the code.

### Files
The R code to run the exposure estimate  
[Exposure_T2_HT2_toGit_131224.Rmd](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Code/Exposure_T2_HT2_toGit_131224.Rmd)

The data file for selection of IDs from the original consumption file  
[Selected_ID.csv](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Selected_ID.csv)

The concentration data obtained from the literature  
[Overview_occurrence_T2HT2.csv](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Overview_occurrence_T2HT2_131224.csv)

The concumption Dummy data  
[Consumption_dummy.csv](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Consumption_dummy.csv)

### Exposure assesment from diet
Establish the folders "Code", "Data" and "Results" in your directory and copy the relevant files to the appropriate folders. Open the R or Rmd code in RStudio and insert your work directory in the code. The data files will be uploaded and create a new folder under "Results" with the date of the day to store the results. 


## Statistical modeling (Rudolf)

A statistical model was developed to describe the relationship between daily intake and excretion of T-2 and HT-2 over time. The model can be interpreted as a compartmental model with transition rates. The input data for the model included urinary concentration data and dietary exposure estimates. The urinary concentration data were assumed to be normally distributed. Youcan read more about the methods in "Exploring the Relationship between Daily Intake and Renal Excretion of the Mycotoxins T-2 and HT-2 Toxin in Humans".  
Individual data from the EuroMix study cannot be openly shared, and we therefore have made dummy data for you to test the code. All the data files below (intake and urinary data) is therfore dummy data. 

### Files

The R code for the statistical modeling  
[Statistical model.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Code/Statistical%20model.R)

Intake data files for the statistical modeling  
[IntakeHT2firstday.xlsx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/IntakeHT2firstday.xlsx)  
[IntakeHT2otherday.xlsx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/IntakeHT2otherday.xlsx)  
[IntakeT2firstday.xlsx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/IntakeT2firstday.xlsx)  
[IntakeT2otherday.xlsx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/IntakeT2otherday.xlsx)  

Urinary data for statistical modelling  
[Urinaryvolumes.csv](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Urinaryvolumes.csv)  
[Urinaryconcentrations.xlsx](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/Urinaryconcentrations.xlsx)  


# Second paper

This is the second paper where we want to study the toxicokinetic of T-2 and HT-2 toxin. In the present study, we wanted to further unravel the toxicokinetics of T-2 and HT-2, aiming to develop a PBK model that can ultimately be used to predict external exposure based on urinary concentrations of the mycotoxins and their metabolites in humans. Eventually, allowing for a more advanced exposure assessment of T-2 and HT-2 and thus, strengthening human risk assessment. This research builds on the work conducted by McKeon et al. (2025) 

## PBPK model developed by J. Westerhout, PhD (joost.westerhout@rivm.nl)

For the development of a PBK model for T-2 and HT-2 we used a population physiology model, namely, physB (Bosgra et al., 2012), as starting point. In short, the physB model provides a set of differential equations describing the kinetics of a chemical in the human body, in which different compartments represent different organs (i.e., stomach, intestine, liver, kidney, lungs, heart, skin, adipose tissue, muscle, brain, gonads, thymus, spleen, pancreas, bone and remaining tissue), which are connected to the systemic circulation by means of arterial and venous blood compartments. Please see the Appendix A for more details. Furthermore, physB includes the variability in physiological characteristics of the human population, which we used to simulate the physiological variability among our participants.

### Files

The R code to run the PBPK model  
[T2_HT2_Phase_I_II_run_model.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Code/T2_HT2_Phase_I_II_run_model.R)  

This is the function files which will be automatically uploaded unto the run file.
[activity.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/functions/activity.R)  
[activityreg.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/functions/activityreg.R)  
[parms_pbk_T2_HT2_Phase_I_II.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/functions/parms_pbk_T2_HT2_Phase_I_II.R)  
[pbk_T2_HT2_Phase_I_II.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/functions/pbk_T2_HT2_Phase_I_II.R)  
[physB.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/functions/physB.R)   
[ptisb_T2_HT2.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/functions/ptisb_T2_HT2.R)  
[utilities.R](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/functions/utilities.R)  

The data that will be used to run the model  
[refman_female.csv](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/refman_female.csv)    
[refman_male.csv](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Data/refman_male.csv)  

The result file with figures  
[output_T2_HT2_Phase_I_II_2026-02-27.Rdata](https://github.com/TrineHusoy/Exposure_T2_HT2/blob/main/Results/output_T2_HT2_Phase_I_II_2026-02-27.Rdata)  


