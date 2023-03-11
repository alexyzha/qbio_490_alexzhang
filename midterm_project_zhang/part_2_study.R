#This will research highest mutation frequency of gene PIK3CA in patients of different age.
##-------------------------------------------------------------------------------------------
#loading packages
if(!require("BiocManager"))  
  install.packages("BiocManager")
if(!require("TCGAbiolinks"))  
  BiocManager::install("TCGAbiolinks")
if(!require("dplyr"))  
  install.packages("dplyr")
if(!require("maftools"))  
  install.packages("maftools")
if(!require("SummarizedExperiment"))  
  install.packages("SummarizedExperiment")
if(!require("ggplot2"))  
  install.packages("ggplot2")
if(!require("DESeq2"))  
  install.packages("DESeq2")
if(!require("survival"))  
  install.packages("survival")
if(!require("tidyr"))  
  install.packages("tidyr")
if(!require("survminer"))  
  install.packages("survminer")
#install.packages(c("data.table", "Biobase", "BiocGenerics", "GenomicRanges"))
if(!require("EnhancedVolcano"))
  install.packages("EnhancedVolcano")
  install.packages("ggrepel")

library(BiocManager)
library(dplyr)
library(ggplot2)
library(TCGAbiolinks)
library(maftools)
library(SummarizedExperiment)
library(DESeq2)
library(survival)
library(tidyr)
library(survminer)
library(data.table)
library(Biobase)
library(BiocGenerics)
library(GenomicRanges)
library(EnhancedVolcano)
library(ggrepel)
##-------------------------------------------------------------------------------------------
#setting wd
setwd('~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/analysis_data')
getwd()
##-------------------------------------------------------------------------------------------
#reading in data
#rna_se has been read in manually from blackboard-downloaded file "rna_se.R"

#reading in clinical.csv
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")

clinic <- GDCprepare_clinic(clin_query,
                            clinical.info = "patient")

clinic <- write.csv(clinic, file='brca_clinical_data.csv', row.names=FALSE)

clinical <- read.csv('brca_clinical_data.csv')

clinical_drug <- GDCprepare_clinic(query = clin_query, clinical.info = "drug") #preparing drug data

clinical_rad <- GDCprepare_clinic(query = clin_query, clinical.info = "radiation") #preparing radiation data

#reading in maf_object
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



#load rna_se related
rna_counts <- rna_se@assays@data$unstranded[,!is.na(rna_se@colData$age_at_index)]
rna_counts <- as.data.frame(rna_counts)

rna_genes <- rna_se@rowRanges@elementMetadata
rna_clinical <- rna_se@colData[!is.na(rna_se@colData$age_at_index), ]
rna_clinical <- as.data.frame(rna_clinical)
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments", F, T)
rna_clinical <- rna_clinical[,treatments_mask]

# name rows and columns of rna_counts using rownames(rna_counts) and colnames(rna_counts)
rownames(rna_counts) <- rownames(rna_counts)
colnames(rna_counts) <- colnames(rna_counts)

tissue_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", F, T)
rna_clinical <- rna_clinical[tissue_mask,]
rna_counts <- rna_counts[, tissue_mask]
row_sums <- rowSums(rna_counts)
low_counts_mask <- ifelse(row_sums > 10, TRUE, FALSE)
rna_counts <- rna_counts[low_counts_mask,]
rna_clinical$age_at_diagnosis <- factor(rna_clinical$age_at_diagnosis)
ageNaMask <- !is.na(rna_clinical$age_at_diagnosis)
vitalStatusMask <- !is.na(rna_clinical$vital_status)

rna_counts <- rna_counts[ ,ageNaMask & vitalStatusMask]
rna_clinical <- rna_clinical[ageNaMask & vitalStatusMask, ]

#load DESeq
?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = rna_counts,
                              colData = rna_clinical,
                              design = ~vital_status + age_category)

?DESeq
dds_obj <- DESeq(dds)]
##-------------------------------------------------------------------------------------------
#WORKING WITH MAF PLOTS
#merge clinical and maf_object@clinical_data
  #clinicalAndMafMerged <- merge(maf_object@clinical.data, clinical, by = "bcr_patient_uuid", multiple = "first")
  #maf_object@clinical.data <- unique(clinicalAndMafMerged)
  #was going to do something with this but my R crashed. After that had trouble loading maf_object and "bcr_patient_uuid" wouldnt even show up on maf_object@clinical.data so idk
  ##update: as of next morning i have fixed the issue
#finding median age to split sample into old and young
median_age <- median(as.numeric(maf_object@clinical.data$age_at_initial_pathologic_diagnosis))
print(median_age) #58 years old
#making category
maf_object@clinical.data$age_category <- ifelse(maf_object@clinical.data$age_at_initial_pathologic_diagnosis > median_age,'Old', 'Young')
#old maf mask
old_mask <- ifelse(maf_object@clinical.data$age_category == 'Old', T, F)
old_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[old_mask]
old_maf <- subsetMaf(maf = maf_object,
                     tsb = old_patient_barcodes)
#young maf mask
young_mask <- ifelse(maf_object@clinical.data$age_category == 'Young', T, F)
young_patient_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[young_mask]
young_maf <- subsetMaf(maf = maf_object,
                     tsb = young_patient_barcodes)
#gene maf
PIK3CA_maf <- subsetMaf(maf = maf_object,
                        genes = 'PIK3CA')

#MAF PLOT ONCOPLOT
oncoplot(maf = maf_object,
         top = 2,
         clinicalFeatures = "age_category",
         borderCol = NA)
ggsave("~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/midterm_project_zhang/outputs/oncoplot_1.png")

#MAF PLOT LOLLIPOP PLOT
lollipopPlot(maf = maf_object,
             gene = "PIK3CA")
ggsave("~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/midterm_project_zhang/outputs/lollipop_plot_1.png")

coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = 'Patients (<= 58)', 
           m2Name = 'Patients (> 58)', 
           borderCol = NA)
ggsave("~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/midterm_project_zhang/outputs/cooncoplot_1.png")

lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name = 'Young Patients (<= 58)',
              m2_name = 'Old Patients (> 58)',
              gene = "PIK3CA")
ggsave("~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/midterm_project_zhang/outputs/lollipop_plot_2.png")
##-------------------------------------------------------------------------------------------
#WORKING WITH KAPLAN MEIER SURVIVAL PLOTS
kpDataMerge <- merge(clinical$vital_status, clinical$age_at_initial_pathologic_diagnosis)
kpMask <- ifelse(kpDataMerge$x == "Alive", 1, 0)
kpDataMerge$x <- as.numeric(kpMask)
#so here, 1 = alive, 0 = dead; will plot relative rates of suvival
kaplan_meier <- survfit(Surv(kpDataMerge$y, kpDataMerge$x) ~1, data = clinical)
#i was dumb and switched x and y around in kpDataMerge
plot(kaplan_meier, xlab = "Age", ylab = "Relative Frequency of Survival", main = "Age and Relative Survival Frequency")
ggsave2("~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/midterm_project_zhang/outputs/kaplan_meier_plot_1.png")
##-------------------------------------------------------------------------------------------
#R PLOT
#similar to Kaplan Meier plots, will be plotting vital status against age
rMerge <- merge(clinical$age_at_initial_pathologic_diagnosis, clinical$vital_status) %>%
  extract(x, into = "Age") %>%
  extract(y, into = "Vitals")
  rMask <- ifelse(clinical$vital_status == "Alive", 1, 0)
  rMerge$vital_status <- rMask
  rMerge$age <- as.numeric(rMerge$Age)
  rMerge$vital <- as.numeric(rMerge$vital_status)
rMerge <- rMerge %>%
  select(-Vitals, -vital_status, -Age)

boxplot(rMerge$age, rMerge$vital, type = "p", main = "Age vs Survival", xlab = "Age", ylab = "Relative Survivability")
ggsave2("~/Desktop/Main/School_Files/S2/QBIO/qbio_490_alexzhang/midterm_project_zhang/outputs/boxplot.png")
##-------------------------------------------------------------------------------------------
#VOLCANO PLOT/DIFFERENTIAL EXPRESSION ANALYSIS
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
PIK3CA_sub <- subset(rna_genes, gene_name == "PIK3CA")
PIK3CA_counts <- as.list(PIK3CA_sub[,-1])

EnhancedVolcano(rna_se, x = "PIK3CA", y = "fold_change", pval = "p_value")

















