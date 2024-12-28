
downloads_path = "C:/Users/Devaharish/Downloads"
file_path = paste(downloads_path, "brca_tcga_pan_can_atlas_2018.tar.gz", sep = "/")
untar(file_path)

folder_path = paste(getwd(), "brca_tcga_pan_can_atlas_2018", sep = "/")

data_patient_path = paste(folder_path, "data_clinical_patient.txt", sep = "/")
data_patient = read.delim(data_patient_path)


data_patient = data_patient[5:dim(data_patient)[1], ]

top_genes <- rownames(res_sig)[1:10]

par(mfrow = c(2, 5)) 
for (gene in top_genes) {
  hist(assay[gene, ], 
       main = paste("Histogram of Expression for", gene), 
       xlab = "Expression Level", 
       col = "lightblue", 
       border = "black")
}

path_RNA = paste(folder_path, "data_mrna_seq_v2_rsem.txt", sep = "/")
data_Rnaseq = read.delim(path_RNA)

assay = round(as.matrix(data_Rnaseq[,-c(1,2)]))
rownames(assay) = data_Rnaseq[,1]

metadata = matrix(0, dim(assay)[2], 2)
pat_ids = data_patient[, 1]
stage_factor = as.factor(stage)

for (i in 1:dim(assay)[2]) {
  pat_barcode = colnames(assay)[i]
  pat_barcode = substr(pat_barcode, 1, 12)
  pat_barcode = gsub("\\.", "-", pat_barcode)
  idx = which(pat_barcode == pat_ids)
  metadata[i, 1] = 1 * (as.numeric(data_patient[idx, col_age]) < 55)  # Early diagnosis (age < 55)
  metadata[i, 2] = stage_factor[idx]
}

metadata[is.na(metadata)] = 0
colnames(metadata) = c("Early", "Stage")

if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
library(DESeq2)

assay[is.na(assay)] = 0 
assay[assay < 0] = 0
smallestGroupSize = 3
keep = rowSums(assay >= 10) >= smallestGroupSize
assay = assay[keep,]

dds = DESeqDataSetFromMatrix(countData = assay, colData = metadata, design = ~ Early + Stage)
dds = DESeq(dds)
resultsNames(dds)

print(DE_genes_data)

DE_genes_filtered <- DE_genes_data[, c("Gene", "log2FoldChange")]

print(DE_genes_filtered)

library(clusterProfiler)
library(org.Hs.eg.db)

DE_over = rownames(res_sig[res_sig$log2FoldChange > 0,])
DE_under = rownames(res_sig[res_sig$log2FoldChange < 0,])

go_results_over <- enrichGO(
  gene = DE_over,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

print(head(go_results_over))

go_results_under <- enrichGO(
  gene = DE_under,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
print(head(go_results_under))

kegg_results_over <- enrichKEGG(
  gene = gene_entrez_over[, 2],  
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

print(head(kegg_results_over))

kegg_results_under <- enrichKEGG(
  gene = gene_entrez_under[, 2],
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

print(head(kegg_results_under))

kegg_results_over <- enrichKEGG(
  gene = gene_entrez_over[, 2], 
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

print(head(kegg_results_over))
kegg_results_under <- enrichKEGG(
  gene = gene_entrez_under[, 2], 
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

reactome_results_over <- enrichPathway(
  gene = gene_entrez_over[, 2], 
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

print(head(reactome_results_over))

reactome_results_under <- enrichPathway(
  gene = gene_entrez_under[, 2],  
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

print(head(reactome_results_under))

dotplot(go_results_over, showCategory = 10) + ggtitle("GO Enrichment Overexpressed Genes")
dotplot(go_results_under, showCategory = 10) + ggtitle("GO Enrichment Underexpressed Genes")
dotplot(kegg_results_over, showCategory = 10) + ggtitle("KEGG Pathway Enrichment Overexpressed Genes")
dotplot(kegg_results_under, showCategory = 10) + ggtitle("KEGG Pathway Enrichment Underexpressed Genes")
dotplot(reactome_results_over, showCategory = 10) + ggtitle("Reactome Pathway Enrichment Overexpressed Genes")
dotplot(reactome_results_under, showCategory = 10) + ggtitle("Reactome Pathway Enrichment Underexpressed Genes")
vsd = vst(dds)

par(mfrow = c(1, 2))
plotPCA(vsd, intgroup = c("Early"))
plotPCA(vsd, intgroup = c("Stage"))

library(pheatmap)

top_DE = order(res$padj)
vsd_DE = assay(vsd)[top_DE[1:20], ]

annotation_colors = list(Early = c(Early = "#1f78b4", Late = "#33a02c"))
annotation_col = data.frame(Early = as.matrix(metadata[, 1]))
rownames(annotation_col) = colnames(vsd)

pheatmap(
  vsd_DE,
  cluster_rows = TRUE,
  cluster_cols = TRUE, 
  annotation_col = annotation_col 
)

library(glmnet)
library(survival)

valid_rows <- !is.na(data_patient$Overall.Survival..Months.) & !is.na(data_patient$Overall.Survival.Status)
data_patient <- data_patient[valid_rows, ]

common_samples <- intersect(data_patient$Patient_ID, colnames(assay(vsd)))
data_patient <- data_patient[data_patient$Patient_ID %in% common_samples, ]
vsd <- vsd[, colnames(assay(vsd)) %in% common_samples]

res_sig <- res_sig[rownames(res_sig) %in% rownames(assay(vsd)), ]

surv_data <- data.frame(
  Survival_Time = as.numeric(data_patient$Overall.Survival..Months.),
  Survival_Status = as.numeric(data_patient$Overall.Survival.Status),
  t(assay(vsd)[rownames(res_sig), ])
)

predictor_matrix <- as.matrix(surv_data[, -c(1, 2)])

predictor_variances <- apply(predictor_matrix, 2, var)
zero_variance_predictors <- which(predictor_variances == 0)

predictor_variances_filtered <- apply(predictor_matrix_filtered, 2, var)
zero_variance_predictors <- which(predictor_variances_filtered == 0)

cat("Number of predictors with zero variance after filtering:", length(zero_variance_predictors), "\n")

non_zero_variance_predictors <- which(predictor_variances_filtered > 0)
cat("Number of predictors with non-zero variance:", length(non_zero_variance_predictors), "\n")

if (length(non_zero_variance_predictors) > 0) {
  predictor_matrix_filtered <- predictor_matrix_filtered[, non_zero_variance_predictors]
  survival_response <- Surv(surv_data$Survival_Time, surv_data$Survival_Status)
  
  cox_model <- cv.glmnet(
    predictor_matrix_filtered, 
    survival_response,
    family = "cox"
  )
  
  print(cox_model)
} else {
  cat("No predictors left with non-zero variance after filtering.\n")
}
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")

library(survival)
library(survminer)

surv_data <- surv_data[!is.na(surv_data$Survival_Time) & !is.na(surv_data$Survival_Status), ]

survival_object <- Surv(surv_data$Survival_Time, surv_data$Survival_Status)

colnames(surv_data)

surv_data$gene_expression <- as.numeric(surv_data$gene_expression)

cox_model <- coxph(survival_object ~ gene_expression, data = surv_data)
summary(cox_model)
surv_fit <- survfit(cox_model)

ggsurvplot(surv_fit, 
           data = surv_data, 
           pval = TRUE, 
           conf.int = TRUE, 
           risk.table = TRUE, 
           ggtheme = theme_minimal(), 
           title = "Survival Analysis by Gene Expression")
