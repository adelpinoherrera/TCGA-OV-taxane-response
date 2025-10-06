# TCGA-OV-taxane-response
Pipeline for analyzing TCGA-OV (ovarian cancer) bulk RNA sequencing data and comparing responses to taxane-based therapies

## Introduction
This repository contains code files (.R) and results files (.csv and ZIP) used to download and analyze TCGA-OV dataset. The code is written in R and can be used in version 4.4 and 4.5.
The data generated included a variety of tables and .rds files that was then input in GraphPad to analyze correlations and comparisons. 

## Rstudio and package requirements 
- Rstudio (4.4 or older)
- dplyr
- readr
- ggplot2
- edgeR
- DESeq2
- tidyr
- TCGAbiolinks
- SummarizedExperiment
- org.Hs.eg.db
- janitor

## Content summary 
- TCGA-OV_taxane-free-intervals_analysis.R code used to download TCGA-OV data both bulk RNAseq and metadata of interest as well as differential gene expression analysis and normalization of transcriptomic information
- TCGA-OV_rawData.rds SummarizedExperiment object that contains both transcriptomic and metadata from the TCGA-OV samples that contained bulk RNAseq information
- allTreatments_219patients_noNA.csv table containing information for all the treatments received by each patients as well as start and end of treatment after deleting rows where one of the previous variables was not available. The table contains a row per treatment received per patient so each patient might contain multiple rows
- 90patients_firstTreatISTax_numOFtax_tfi.csv table containing 1 per patient with the following information: arcode, days_to_treatment_start, days_to_treatment_end, therapeutic_agents, taxane_first, taxane_second, number_taxanes, taxane_type, tfi
- TCGA-OVCA_count_matrix.zip table with the bulkRNAseq count matrix for all the patients in the TCGA-OV dataset with gene symbols
- differetially_expressed_genes-HigherYESVSLowerNO-secondTax_90patients.csv table that contains GeneSymbol, log2FoldChange, pvalue and padj for the differential expressed genes between the patients that received a taxane-based therapeutic as a second treatment compared to those that received an alternative treatment
-  processed_data_90patients_TFI-from-firstTaxaneANDnextTreat_logNorm.zip log normalized matrix for the 90 patients
-  lognorm_90patients_count_matrix_WITH_tfi_firstTaxaneANDnexTreat.zip transposed log normalized matrix for the 90 patients with metadata
  
