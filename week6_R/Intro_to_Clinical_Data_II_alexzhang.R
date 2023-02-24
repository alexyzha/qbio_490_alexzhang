##------------------------------------------------------------------------------------------------
#loading packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install(version = "3.16") #install BiocManager
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks") #install TCGAbiolinks
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools") #install maftools
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr") #install dplyr
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr") #install tidyr
if (!require("survival", quietly = TRUE))
  install.packages("survival") #install survival
  
library(BiocManager)
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(tidyr)
library(survival)
#library to load all packages
##------------------------------------------------------------------------------------------------
#setting directory
setwd('/Users/aly/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/analysis_data') #set directory
getwd() #double check directory
##------------------------------------------------------------------------------------------------
#answering the following questions as comments
clinic_data_csv <- read.csv('brca_clinical_data.csv') #read csv
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")
#clin_query
clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")
#clinic prepare
clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug") #preparing drug data
clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation") #preparing radiation data

#answering the following questions
colnames(clinic)
#1. chose "days_to_birth"
is.na(clinic$days_to_birth) #not too many NA
data.class(clinic$days_to_birth)
#2. it is numeric data class therefore continuous
colnames(clinical_drug) #looking through variables for clinical_drug
colnames(clinical_rad) #looking thorugh variables for clinical_rad
data.class(clinical_drug$total_dose) #found one, looking at data class
#4. it says factor, therefore categorical, as factor separates into numeric or strings
clinical_drug$total_dose #looking at the variable
clinical_drug$total_dose_units #loaded in related as dose itself isn't standardized for units. ex: 1 ug and 1 mg both = 1.
#3. I chose the variable total_dose from clinical_drug. there is a huge chunk of data with NAs, the data is numeric as far as i can see
#5. 1. idk if it is possible to plot total_dose against days_to_birth, as the columns from the 2 dataframes have different numbers of obs.
# hypothesis 1: as age (days_to_birth) increases, total_dose decreases
# hypothesis 2: as age increases, survival in breast cancer decreases
# hypothesis 3: as total_dose increases, survival in breast cancer increases
##------------------------------------------------------------------------------------------------
#coding activity
#1
  #ignore this chunk of comments
    #ageMask <- !is.na(clinic$days_to_birth) #testing masking
    #ageMASKED <- clinic$days_to_birth[ageMask] 
    #doseMask <- !is.na(clinical_drug$total_dose)
    #doseMASKED <- clinical_drug$total_dose[doseMask]
    #plot(clinic$days_to_birth, clinical_drug$total_dose, type = "p", main = "test", xlim = c(10000, 30000), ylim = c(0, 3000)) #test, can't compare the 2 variables bc different lengths
  #package installation (ignore also)
    #install.packages("dplyr") #install dplyr
    #install.packages("tidyr") #install tidyr
    #library(dplyr) #loading in dplyr
    #library(tidyr) #loading tidyr
#making plot below
mergedForPlot <- inner_join(clinic, clinical_drug, by = "bcr_patient_barcode", multiple = "first") %>%
  dplyr::select(days_to_birth, total_dose) %>%
  extract(days_to_birth, into = "numAGE") %>%
  extract(total_dose, into = "numDOSE") %>%
  dplyr::mutate(as.numeric(numAGE)) %>%
  dplyr::mutate(as.numeric(numDOSE))
#whole thing above should be run together bc pipe
plot(mergedForPlot$numAGE, mergedForPlot$numDOSE, type = "p", main = "Age vs Total Dose", xlim = c(10000, 30000), ylim = c(0, 3000), xlab = "Age (days)", ylab = "Total Dose") #test, comparing with merged dataframes
#plotted the 2 against each other. one issue is that I couldn't figure out how to standardize all dosage values by mg; there are some with ug measurements which may be outliers
plot(mergedForPlot$numAGE, mergedForPlot$numDOSE, type = "h", main = "Age vs Total Dose", xlim = c(10000, 30000), ylim = c(0, 3000), xlab = "Age (days)", ylab = "Total Dose") #test, comparing with merged dataframes
#histogram above
#from the plot functions above, its evident that the dosage amount is generally very low regardless of age. However, with higher age, it seems to decrease slightly. ~20000 days of age is when dosage is highest, but it may also be where the most common ages are

#2
#install.packages("survival") #install survival pkg
#library(survival) #load survival for Kaplan-Meier
#looking at vital_status and days_to_birth to see what i have to do to it
  clinic$vital_status
  data.class(clinic$vital_status)
#bc ill need the merge later for alive/dead status, i merge again
  mergedForCode <- inner_join(clinic, clinical_drug, by = "bcr_patient_barcode", multiple = "first")
  data.class(mergedForCode$days_to_birth) #can only see mergedForCode$days_to_birth after defining mergedForCode
#now i convert vital_status into 1's and 0's, adding 1 column to mergedForCode:
  mergedForCode$numVS <- ifelse(mergedForCode$vital_status == "Alive", 1, 0) #to make alives 1 and deads 0
  mergedForCode$numVS #to check if it changed (it did)
  mergedForCode$days_to_birth <- ifelse(mergedForCode$days_to_birth < 0, mergedForCode$days_to_birth * -1, mergedForCode$days_to_birth) #to make negatives positive
  mergedForCode$days_to_birth #to check if it changed (it did)
#making Kaplan-Meier plot 1
  KPfit1 <- survfit(Surv(mergedForCode$days_to_birth, mergedForCode$vital_status) ~1, data = mergedForCode)
  plot(KPfit1, xlab = "Age (days)", ylab = "Relative Frequency of Survival", main = "Age and Relative Survival Frequency")

#3
#packages installed, looking at data
  mergedForCode
  mergedForCode$total_dose #looking at mergedForCode$total_dose
  mergedForCode$total_dose <- as.numeric(mergedForPlot$numDOSE) #making everything numeric
  data.class(mergedForCode$total_dose) #checking to make sure its numeric (it is)
  mergedForCode$total_dose #check to make sure it is correct
#making Kaplan-Meier plot 2
  KPfit2 <- survfit(Surv(mergedForCode$total_dose, mergedForCode$vital_status) ~1, data = mergedForCode)
  plot(KPfit2, xlab = "Total Dosage (mg)", ylab = "Relative Frequency of Survival", main = "Total Dosage and Relative Survival Frequency")
  
  
  
  
  
  
  
  
  
  
  
  

