path <- "C:/Users/Devaharish/Downloads"
folder_name <- "brca_tcga_pan_can_atlas_2018.tar.gz"
folder <- paste(path, folder_name, sep = "/")

untar(folder)

new_dir <- paste(getwd(), "brca_tcga_pan_can_atlas_2018", sep = "/")
setwd(new_dir)

data_patient <- read.delim("data_clinical_patient.txt")  
data_Rnaseq <- read.delim("data_mrna_seq_v2_rsem.txt")  
data_CNA <- read.delim("data_cna.txt")                 

assay <- as.matrix(data_Rnaseq[,-c(1,2)]) 

metadata <- data.frame(ERBB2_Amplified = rep(0, ncol(assay)))

pat_ids <- data_patient[, 1]

for (i in 1:ncol(assay)) {
  pat_barcode <- colnames(assay)[i]
  pat_barcode <- substr(pat_barcode, 1, 12) 
  pat_barcode <- gsub("\\.", "-", pat_barcode) 
  
  idx_cna <- which(data_CNA$SAMPLE_ID == pat_barcode)
  if (length(idx_cna) > 0 && data_CNA$ERBB2[idx_cna] > 0) {
    metadata$ERBB2_Amplified[i] <- 1  
  }
}

metadata$ERBB2_Amplified[is.na(metadata$ERBB2_Amplified)] <- 0

metadata$ERBB2_Amplified <- as.factor(metadata$ERBB2_Amplified)

str(metadata)

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")
library(DESeq2)

assay[is.na(assay)] <- 0
assay[assay < 0] <- 0 

dds <- DESeqDataSetFromMatrix(countData = round(assay),
                              colData = metadata,
                              design = ~ ERBB2_Amplified)

dds <- DESeq(dds)

vsd <- vst(dds)
head(assay(vsd))

