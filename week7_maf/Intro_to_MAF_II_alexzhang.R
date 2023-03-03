##------------------------------------------------------------------------------------------
#set directory
setwd('~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/analysis_data')
#check
getwd()
##------------------------------------------------------------------------------------------
#package loading
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16") #install BiocManager
if (!require("maftools", quietly = TRUE))
  BiocManager::install("maftools") #install maftools
if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks") #install TCGAbiolinks
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2") #install ggplot2

library(maftools)
library(ggplot2)
library(TCGAbiolinks)
library(BiocManager)
##------------------------------------------------------------------------------------------
#loading in clinical data
clinical <- read.csv("/Users/aly/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/analysis_data/brca_clinical_data.csv")

#loading maf_object
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(maf_query)

maf <- GDCprepare(maf_query)

maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)
##------------------------------------------------------------------------------------------
#questions 1-5
#q1, using gender
maf_object@clinical.data$gender_category <- ifelse(maf_object@clinical.data$gender == "MALE", "MALE", "FEMALE")
#male mask
male_mask <- ifelse(maf_object@clinical.data$gender_category == "MALE", T, F)
male_patient_barcodes <- maf_object@clinical.data$gender_category[male_mask]
#female mask
female_mask <- ifelse(maf_object@clinical.data$gender_category == "FEMALE", T, F)
female_patient_barcodes <- maf_object@clinical.data$gender_category[female_mask]

#q2
oncoplot(maf = maf_object,
         top = 10,
         clinicalFeatures = "gender_category",
         borderCol = NA)
#save
ggsave("/Users/aly/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/week7_maf/part_2_oncoplot.png")
#top 10 mutated genes for the 2 groups are TP53, PIK3CA, TTN, CDH1, GATA3, MUC16, KMT2C, MAP3K1, HMCN1, and FLG
#going to use PIK3CA
#PIK3CA is a gene plays a role in signaling pathways that regulate cell growth, proliferation, and survival.
#PIK3CA mutations were found to be more common in women than in men, as they are associated with an increased risk of developing breast, ovarian and endometrial cancers. Also, PIK3CA mutations may be associated with an increased risk of recurrence in breast cancer patients.

#a3
#define vars for contingency
gender_maf <- subsetMaf(maf = maf_object,
                     tsb = maf_object@clinical.data$gender_category)
PIK3CA_maf <- subsetMaf(maf = maf_object,
                        gene = 'PIK3CA')
#make contingency table
contig <- table(gender_maf, 
                PIK3CA_maf)
#plot contingency table
mosaicplot(contig)
#save
ggsave("/Users/aly/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/week7_maf/part_2_contingency.png")

#q4
#define male maf
male_maf <- subsetMaf(maf = maf_object,
                      tsb = male_patient_barcodes)
#define female maf
female_maf <- subsetMaf(maf = maf_object,
                        tsb = female_patient_barcodes)
#making cooncoplot
coOncoplot(m1 = male_maf, 
           m2 = female_maf, 
           m1Name = 'Male Patients', 
           m2Name = 'Female Patients', 
           borderCol = NA)
#save
ggsave("/Users/aly/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/week7_maf/part_2_cooncoplot")
#the issue i found with using gender was that there wasn't enough data from the male samples to tell anything really. There were only 11 male samples...

#q5
maf_object@clinical.data$Overall_Survival_Status <- ifelse(maf_object@clinical.data$gender == "MALE", T, F)

mafSurvival(maf = maf_object,
            genes = "PIK3CA",
            time = "days_to_last_followup",
            Status = "Overall_Survival_Status",
            isTCGA = TRUE)
#error in plot.new(): figure margins too large
#i think there isn't enough male data to do this with :(
#i started the assignment too late to switch variables from gender but i think i have a pretty good understanding on the code needed to do oncoplots, contingency plots, cooncoplots, and survival curves
#pls be nice to me :( it was midterm week















