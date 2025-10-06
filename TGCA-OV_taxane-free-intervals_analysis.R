#entire workflow for TCGA-taxane responses
###Load packages
library('dplyr') 
library('readr')
library(ggplot2)
library(edgeR)
library(DESeq2)
library(tidyr)

BiocManager::install("TCGAbiolinks") #specific to download TCGA data 

library('TCGAbiolinks') #only for downloading data 
library('SummarizedExperiment') #only for downloading data 

library(org.Hs.eg.db) #packages for converting ENSEMBL symbols to gene symbols
library(janitor)

#### Download TCGA-OV dataset, data downloaded and saved it ####
#only downloaded the transcriptomic data and includes all the clinical data
OVCA_query <- GDCquery(project = "TCGA-OV", #specofy name of the dataset wanted 
                       data.category = "Transcriptome Profiling", #specify bulkRNAseq data 
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts")
#download dataset 
GDCdownload(OVCA_query)
#get results for the dataset 
OVCA_res = getResults(OVCA_query)

#make a summarizedExperiment object and save
OVCA_data <- GDCprepare(OVCA_query)
saveRDS(OVCA_data, "/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/TCGA-OV/taxane-free-intervals/publication/TCGA-OV_rawData.rds")


####Create a dataframe with information of interest ###
#extract all the metadata
col_data <- colData(OVCA_data)
#extract patients barcodes
barcodes <- col_data$barcode #429 barcodes for patients 
#extract days to diagnosis
days_to_diagnosis <- col_data$days_to_diagnosis
#extract all the treatments that the patients received, start and end of treatment dates
treatment_list <- col_data$treatments
#prior_treatment <- col_data$prior_treatment #originally 278 NO and 149 YES 
#days_to_last_follow_up <- col_data$days_to_last_follow_up

###Now unpack treatments to get a row per treatment received (multiple treatments per patients) start and end of treatment
all_treatments <- lapply(seq_along(treatment_list), function(i) {
  treatments <- treatment_list[[i]]
  
  #extract treatments and skip if no treatment was received
  if (is.null(treatments) || nrow(treatments) == 0) return(NULL)
  
  n <- nrow(treatments)
  
  #extract days to treatment start and skip if no treatment start 
  days_start <- treatments$days_to_treatment_start
  if (is.null(days_start) || length(days_start) != n) {
    days_start <- rep(NA, n)
  }
  
  #extract days to treatment end and skip if no treatment end 
  days_end <- treatments$days_to_treatment_end
  if (is.null(days_end) || length(days_end) != n) {
    days_end <- rep(NA, n)
  }
  
  #extract each therapeutic agent received
  agents <- treatments$therapeutic_agents
  if (is.null(agents) || length(agents) != n) {
    agents_collapsed <- rep(NA_character_, n)
  } else if (is.list(agents)) {
    agents_collapsed <- sapply(agents, function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      paste(x, collapse = ", ")
    })
  } else {
    agents_collapsed <- as.character(agents)
  }
  
  # Create the data frame with matching row counts
  data.frame(
    barcode = rep(barcodes[i], n),
    days_to_diagnosis = rep(days_to_diagnosis[i], n),
    days_to_treatment_start = days_start,
    days_to_treatment_end = days_end,
    therapeutic_agents = agents_collapsed,
    stringsAsFactors = FALSE
  )
})

##create table with one row per therapeutic agent received by the patients including days_to_treatment_start and days_to_treatment_end
  ##the table should include multiple rows for the same patient to include each therapeutic agent they received 
treatment_df2 <- do.call(rbind, all_treatments)
unique(treatment_df2$barcode) #406 patients

###get the table with the patients that have all treatment information variable
##delete rows with NA entries in any of the following columns: days_to_treatment_start, days_to_treatment_end, therapeutic_agents
treatment_df_filtered <- treatment_df2 %>%
  filter(
    !is.na(days_to_treatment_start) &
      !is.na(days_to_treatment_end) &
      !is.na(therapeutic_agents)
  )

unique(treatment_df_filtered$barcode) #219 patients with treatment information 

#save the table
write_csv(treatment_df_filtered, "/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/TCGA-OV/taxane-free-intervals/publication/allTreatments_219patients_noNA.csv")
#open table on excel keep patients that received taxane as their first treatment
  #calculate taxane-free interval (TFI) from days of treatment end of the first treatment to the days of treatment start of the second treatment 
    #find number of taxanes received per patient 
      #note if the second treatment was taxane-based or not 
        #make the final table have 1 per patient with the following entries: barcode, days_to_treatment_start, days_to_treatment_end, therapeutic_agents
        #taxane_first, taxane_second, number_taxanes, taxane_type, tfi

#### Upload the table again with 90 patients #### 
patient90 <- read.csv("/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/TCGA-OV/taxane-free-intervals/90patients_firstTreatISTax_numOFtax_tfi.csv")

##create new category columns >180 and <180 and >90 and <90 
patient90 <- patient90 %>%
  mutate(tfi_180group = ifelse(tfi > 180, "above_180", "below_180"))
table(patient90$tfi_180group)
# above_180 below_180 
# 40        50

patient90 <- patient90 %>%
  mutate(numberTax_group = ifelse(number_taxanes >= 2, "aboveORequal_2", "below_2"))
table(patient90$numberTax_group)
#aboveORequal_2 below_2 
#49      41 


#### Get the count matrix for the 90 patients ####
##obtain the count matrix with counts from the summarized experiment object
OVCA_count_matrix <- assay(OVCA_data)
##transform to a dataframe 
#OVCA_count_matrix_df <- as.data.frame(OVCA_count_matrix)
  #genes on the count matrix in the ENSEMBl symbols so change to gene symbols

##remove the version number at the end of the ENSEMBL id 
OVCA_genes <- rownames(OVCA_count_matrix_df) %>%
  tibble::enframe() %>%
  mutate(ENSEMBL = stringr::str_replace(value, "\\.[0-9]+",""))


##change from ENSEMBL to symbol 
#keep only one gene out of the duplicates
OVCA_genes_map <- clusterProfiler::bitr(OVCA_genes$ENSEMBL,
                                        fromType = "ENSEMBL",
                                        toType = "SYMBOL",
                                        OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, .keep_all = TRUE)
#all the symbols in the same dataframe to double check on any errors
OVCA_genes_map <- OVCA_genes_map %>%
  left_join(OVCA_genes)

#reduce count matrix to only include the ENSEMBL symbols that were able to be converted
OVCA_count_matrix <- OVCA_count_matrix[OVCA_genes_map$ENSEMBL, ]
#change rownames to include the gene symbols 
row.names(OVCA_count_matrix) <- OVCA_genes_map$SYMBOL
#convert to a dataframe 
OVCA_count_matrix_df <- data.frame(OVCA_count_matrix)
#save the file with all the patients 
write.csv(OVCA_count_matrix_df,'/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/TCGA-OV/taxane-free-intervals/publication/TCGA-OVCA_count_matrix.csv',row.names=TRUE)

##reduce count matrix to only include the 90 patients of interest 
#change the colnames to have - between the numbers and not .
colnames(OVCA_count_matrix_df) <- gsub("\\.", "-", colnames(OVCA_count_matrix_df))

#collect barcodes for the 90 patients 
barcodes_90patients <- patient90$barcode
#reduce the count matrix to only contain the information for the 90 patients of interest 
count_matrix_reduced <- OVCA_count_matrix_df[,colnames(OVCA_count_matrix_df) %in% barcodes_90patients] #90 columns 
dim(count_matrix_reduced) #37464   90


#### Normalize this new matrix and do differential gene expression #### 
#calculate standard deviation across all samples, no need to exclude the first column because the gene symbols are the rownames
count_matrix_sd <- count_matrix_reduced[order(-apply(count_matrix_reduced, 1, sd)), ]
dim(count_matrix_sd) #37464   90

#remove duplicated genes, in this case all of this was done prior 
count_matrix_sd <- count_matrix_sd[!duplicated(rownames(count_matrix_sd)) ,] #37464 leftover genes, no duplicated genes
#remove genes that were not mapped 
count_matrix_sd <- count_matrix_sd[!is.na(rownames(count_matrix_sd)),] #37464 leftover genes, all genes were mapped


#convert data (excluding the first column) to a numeric matrix 
count_matrix_sd1 <- as.matrix(count_matrix_sd)

#filtering the count matrix, based on a threshold for CPM
minCPM <- 0.5 
nLibraries <- 1
count_matrix_sd_matrix_filtered <- count_matrix_sd1[ which( apply( cpm(DGEList(counts = count_matrix_sd1)),  1,
                                                                  function(y) sum(y >= minCPM) ) >= nLibraries ), ] 
dim(count_matrix_sd_matrix_filtered) #22721 90 after filtering 

###build DESeq2, to compare data in the second taxane group 
##DESeq2 needs colData to be able to run

#make a metadata factor based on tfi_group
rownames(patient90) <- patient90$barcode
#check that all the patients are both in the 90patient (metadata) count_matrix_sd_matrix_filtered (count matrix) objects 
all(colnames(count_matrix_sd_matrix_filtered) %in% rownames(patient90)) #true 
#check that they are in the same order if not it won't run 
all(rownames(patient90) == colnames(count_matrix_sd_matrix_filtered)) #false, not in the same order 
#re order columns in count matrix to match the metadata matrix
patient90 <- patient90[colnames(count_matrix_sd_matrix_filtered), ]
all(rownames(patient90) == colnames(count_matrix_sd_matrix_filtered)) #true, in the same order 

#run dds based on receiving a taxane-based treatment as the second treatment or not 
  #design variable needs to be present in the 90patient object 
dds <- DESeqDataSetFromMatrix(countData = count_matrix_sd_matrix_filtered, colData = patient90, design = ~ taxane_second) 
dds <- DESeq(dds)
res <- results(dds, contrast = c("taxane_second", "yes","no"))
#create table with gene names, log2FC, pval and padj 
res$GeneSymbol <- rownames(res)
#reorder columns: 
res_df <- res[,c("GeneSymbol", "log2FoldChange", "pvalue", "padj")]
res_df_table <- as.data.frame(res_df@listData)
#save the results 
write.csv(res_df_table,"/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/TCGA-OV/taxane-free-intervals/publication/differetially_expressed_genes-HigherYESVSLowerNO-secondTax_90patients.csv" )


#lastly do log transformation 
c <- 4 #pseudocount is 4
count_matrix_log_norm <- log2(counts(dds, normalized=TRUE) + c)
write.csv(count_matrix_log_norm,"/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/TCGA-OV/taxane-free-intervals/publication/processed_data_90patients_TFI-from-firstTaxaneANDnextTreat_logNorm.csv")

#### Attempt to add metadata to the log normalized count matrix table ####
#transpose matrix
count_matrix_log_norm_t <- t(count_matrix_log_norm)

#add a barcode column on the transpose count matrix 
count_matrix_log_norm_t <- data.frame(barcode = rownames(count_matrix_log_norm_t), count_matrix_log_norm_t)

#delete non-relevant information from 90 patient dataframe
patient90$days_to_treatment_start <- NULL
patient90$days_to_treatment_end <- NULL
patient90$therapeutic_agents <- NULL

#merge count matrix with metadata table
count_matrix_with_90patients <- merge(patient90, count_matrix_log_norm_t, by = "barcode", all.x = TRUE)
write_csv(count_matrix_with_90patients, '/blue/ferrallm/01_analysis/adelpinoherrera-bulkRNAseq/TCGA-OV/taxane-free-intervals/publication/lognorm_90patients_count_matrix_WITH_tfi_firstTaxaneANDnexTreat.csv')

###use the values for TFI and expression for markers of interest (CFLAR, ABCB1 and STAT3) to create correlations using graphpad 


